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

  // Filter Error Messages
  void PrintEKFProgressHeaderPVC() {
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time [min]  | Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms]" << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;
}

VecXd ComputeEstimationErrorPVC(const Ptr<Spacecraft> sat, IFilter* filter) {
  auto x_est = filter->GetUpdatedStateEstimate();
  auto x_true = sat->GetStateVec();

  double x_pos_err = 1000 * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double x_vel_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double x_clk_bias_err = 3e8 * abs((x_true(6) - x_est(6)).val());
  double x_clk_drift_err = 3e8 * abs((x_true(7) - x_est(7)).val());

  VecXd est_err(4);
  est_err << x_pos_err, x_vel_err, x_clk_bias_err, x_clk_drift_err;

  return est_err;
}

void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err) {
  std::cout.precision(5);
  std::cout << std::left << std::setw(12) << t / 60 << " " << std::left << std::setw(12)
            << x_pos_err << "  " << std::left << std::setw(14) << x_vel_err << "   " << std::left
            << std::setw(16) << x_clk_bias_err << std::endl;
};






}  // namespace lupnt
