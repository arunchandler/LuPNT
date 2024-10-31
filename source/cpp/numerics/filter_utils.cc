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

  MatXd ConstructInitCovarianceRVC(double pos_err, double vel_err, double clk_bias_err,
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

  FilterProcessNoiseFunction ConstructProcessNoiseRVC(ClockModel cmodel, int state_size,
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

}  // namespace lupnt
