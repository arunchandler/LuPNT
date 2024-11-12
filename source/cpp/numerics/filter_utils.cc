/**
 * @file filter_utils.cc
 * @author Stanford NAV Lab
 * @brief Utility functions for filters
 * @version 0.1
 * @date 2024-10-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/numerics/filter_utils.h"

namespace lupnt {

  MatXd ConstructInitCovariancePVC(double pos_err, double vel_err, double clk_bias_err,
                                   double clk_drift_err) {
    Mat6d P_rv = Mat6d::Zero();
    P_rv.block(0, 0, 3, 3) = Mat3d::Identity() * pow(pos_err, 2);
    P_rv.block(3, 3, 3, 3) = Mat3d::Identity() * pow(vel_err, 2);

    Mat2d P_clk = Mat2d::Zero();
    P_clk(0, 0) = pow(clk_bias_err, 2);
    P_clk(1, 1) = pow(clk_drift_err, 2);

    MatXd P0 = BlkDiagD(P_rv, P_clk);

    return P0;
  };

  FilterProcessNoiseFunction ConstructProcessNoisePVC(ClockModel cmodel, int state_size,
                                                      double sigma_acc) {
    FilterProcessNoiseFunction proc_noise_func
        = [cmodel, state_size, sigma_acc](const VecX x, Real t_curr, Real t_end) -> MatXd {
      int clock_index = 6;
      double dt = (t_end - t_curr).val();

      MatXd Q = MatXd::Zero(state_size, state_size);

      Mat6d Q_rv = Mat6d::Zero();
      for (int i = 0; i < 3; i++) {
        Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
        Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
        Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
        Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
      }

      Mat2d Q_clk = ClockDynamics::TwoStateNoise(cmodel, dt).cast<double>();

      Q.block(0, 0, 6, 6) = Q_rv;
      Q.block(6, 6, 2, 2) = Q_clk;

      return Q;
    };

    return proc_noise_func;
  }

  FilterProcessNoiseFunction ConstructProcessNoisePVC(ClockModel cmodel, int state_size, double sigma_acc, int n_sat){
    FilterProcessNoiseFunction proc_noise_func
        = [cmodel, state_size, sigma_acc, n_sat](const VecX x, Real t_curr, Real t_end) -> MatXd {
      int clock_index = 6;
      double dt = (t_end - t_curr).val();

      MatXd Q = MatXd::Zero(state_size*n_sat, state_size*n_sat);

      for (int k = 0; k < n_sat; k++){
        Mat6d Q_rv = Mat6d::Zero();
        for (int i = 0; i < 3; i++) {
          Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
          Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
          Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
          Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
        }

        Mat2d Q_clk = ClockDynamics::TwoStateNoise(cmodel, dt).cast<double>();

        Q.block(k*state_size, k*state_size, 6, 6) = Q_rv;
        Q.block(k*state_size+6, k*state_size+6, 2, 2) = Q_clk;
      }

      return Q;
    };

    return proc_noise_func;
  }

// Filter Error Messages
void PrintEKFProgressHeaderPVC() {
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time [min]  | Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms]" << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;
}

void PrintEKFProgressHeaderPVC(int n_sat){
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  
  //        |                  Sat1                     |                   Sat2                    | Sat3  ...
  //  Time  | PosErr [m] | VelErr [mm/s] | ClkBErr [ms] | PosErr [m] | VelErr [mm/s] | ClkBErr [ms] | ...
  // -------------------------------------------------------------------------------------------------------

  std::cout << "            |";
  // first row, plot each sat
  for (int i = 0; i < n_sat; i++){
    std::cout << "Sat" << i+1 << "                                             | ";
  }
  std::cout << std::endl;

  // second row, time and pos, vel and clk bias error for each sat
  std::cout << "Time [min]  |";
  for (int i = 0; i < n_sat; i++){
    std::cout << "Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms] | ";
  }
  std::cout << std::endl;

  for (int i = 0; i < n_sat; i++){
    std::cout << "----------------------------------------------------------";
  }
  std::cout << std::endl;

}

VecXd ComputeEstimationErrorPVC(const Ptr<Spacecraft> sat, IFilter* filter, int start_idx) {
  auto x_est = filter->GetUpdatedStateEstimate();
  auto x_true = sat->GetStateVec();


  double x_pos_err = 1000 * (x_true.segment(0, 3) - x_est.segment(start_idx, 3)).norm().val();
  double x_vel_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(start_idx+3, 3)).norm().val();
  double x_clk_bias_err = 3e8 * abs((x_true(6) - x_est(start_idx+6)).val());
  double x_clk_drift_err = 3e8 * abs((x_true(7) - x_est(start_idx+7)).val());

  VecXd est_err(4);
  est_err << x_pos_err, x_vel_err, x_clk_bias_err, x_clk_drift_err;

  return est_err;
}

VecXd ComputeEstimationErrorPVC(const std::vector<Ptr<Spacecraft>>& sats, IFilter* filter) {

  VecXd est_err = VecXd::Zero(4 * sats.size());

  for (int i = 0; i < sats.size(); i++) {
    VecXd est_err_i = ComputeEstimationErrorPVC(sats[i], filter, i * 8);
    est_err.segment(i * 4, 4) = est_err_i;
  }

  return est_err;
}

VecX ConstructTrueStateVecFromSats(const std::vector<Ptr<Spacecraft>>& sats) {
  VecX x_true = VecX::Zero(8 * sats.size());
  for (int i = 0; i < sats.size(); i++) {
    VecX x_sat = sats[i]->GetStateVec();
    x_true.segment(i * 8, 8) = x_sat;
  }

  return x_true;
}

void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err) {
  std::cout.precision(5);
  std::cout << std::left << std::setw(12) << t / 60 << " " << std::left << std::setw(12)
            << x_pos_err << "  " << std::left << std::setw(14) << x_vel_err << "   " << std::left
            << std::setw(16) << x_clk_bias_err << std::endl;
};

void PrintEKFProgressPVC(double t, const VecXd& est_err, int n_sat) {
  std::cout.precision(5);
  std::cout << std::left << std::setw(12) << t / 60 << " ";
  // for each satellite
  for (int i = 0; i < n_sat; i++){
    std::cout << std::left << std::setw(12) << est_err(i*4) << "  " << std::left << std::setw(14) << est_err(i*4+1) << "   " << std::left << std::setw(16) << est_err(i*4+2) << "  | ";
  }
  std::cout << std::endl;

};

void PrintEstimationStatistics(const VecXd& num_meas, const MatXd& error_mat, double ratio, int n_sat) {
  int n_time = num_meas.size();
  Vec4d rms, means, stds, p68, p95, p99;

  if (error_mat.rows() != 4 * n_sat) {
    std::cout << "Wrong Mat Size, Error Mat size must be (4 x timestep)" << std::endl;
    return;
  }

  // extract statistics range data
  int start_idx = (int)((1.0 - ratio) * n_time);
  int end_idx = n_time - 1;
  int n_range = end_idx - start_idx;

  VecXd num_meas_range(n_range);
  MatXd error_mat_range(4 * n_sat, n_range);

  num_meas_range = num_meas.segment(start_idx, n_range);
  error_mat_range = error_mat.block(0, start_idx, 4 * n_sat, n_range);

  // reshape error_mat_range from (4xnsat, n_range) to (4, n_sat*n_range)
  MatXd error_mat_range_reshaped(4, n_sat * n_range);
  for (int i = 0; i < n_sat; i++) {
    error_mat_range_reshaped.block(0, i * n_range, 4, n_range) = error_mat_range.block(i * 4, 0, 4, n_range);
  }

  // compute statistics ----------------------------------------------------
  // rms
  for (int i = 0; i < 4; i++) {
    rms(i) = RootMeanSquareD(error_mat_range_reshaped.row(i));
    means(i) = error_mat_range_reshaped.row(i).mean();
    stds(i) = StdD(error_mat_range_reshaped.row(i));
    p68(i) = PercentileD(error_mat_range_reshaped.row(i), 0.68);
    p95(i) = PercentileD(error_mat_range_reshaped.row(i), 0.95);
    p99(i) = PercentileD(error_mat_range_reshaped.row(i), 0.99);
  }

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "< Simulation Statistics (Last " << ratio * 100 << "%)>" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Statistics  | Position [m]  | Velocity [mm/s] | Clock Bias [ns] "
               "| Clk Drift [ns/s] "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------------"
            << std::endl;

  std::cout.precision(5);
  std::cout << "RMS         | " << std::left << std::setw(16) << rms(0) << "  " << std::left
            << std::setw(16) << rms(1) << "   " << std::left << std::setw(16) << rms(2) << std::left
            << std::setw(16) << rms(3) << std::endl;

  std::cout << "Mean+-Std   | " << std::left << means(0) << "+-" << std::left << stds(0) << "  "
            << std::left << means(1) << "+-" << std::left << stds(1) << "   " << std::left
            << means(2) << "+-" << std::left << stds(2) << "   " << std::left << means(3) << "+-"
            << std::left << stds(3) << std::endl;

  std::cout << "68%         | " << std::left << std::setw(16) << p68(0) << "  " << std::left
            << std::setw(16) << p68(1) << "   " << std::left << std::setw(16) << p68(2) << std::left
            << std::setw(16) << p68(3) << std::endl;

  std::cout << "95%         | " << std::left << std::setw(16) << p95(0) << "  " << std::left
            << std::setw(16) << p95(1) << "   " << std::left << std::setw(16) << p95(2) << std::left
            << std::setw(16) << p95(3) << std::endl;

  std::cout << "99%         | " << std::left << std::setw(16) << p99(0) << "  " << std::left
            << std::setw(16) << p99(1) << "   " << std::left << std::setw(16) << p99(2) << std::left
            << std::setw(16) << p99(3) << std::endl;

  std::cout << "  " << std::endl;

}






}  // namespace lupnt
