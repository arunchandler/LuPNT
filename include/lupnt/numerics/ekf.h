/**
 * @file filters.h
 * @author Stanford NAV LAB
 * @brief List of Filters
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <Eigen/QR>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/filters.h"

namespace lupnt {
  class EKF : public IFilter {
  public:
    Real t_;
    VecX x_;
    MatXd P_;

    VecX x_prior_;
    VecX x_post_;
    MatXd P_prior_;
    MatXd P_post_;

    MatXd F_;
    MatXd H_;
    MatXd Q_;  // Process noise cov
    MatXd R_;  // Measurement noise cov

    VecX dy_;       // Measurement residual
    VecX dx_;       // State update
    VecX z_true_;   // Observed measurement
    VecX z_prior_;  // Predicted measurement

    MatXd S_;  // Innovation cov
    MatXd K_;  // Kalman gain

    FilterDynamicsFunction f_dyn_;
    FilterProcessNoiseFunction f_proc_;
    FilterMeasurementFunction f_meas_;

    double outlier_threshold_ = 3.0;

    void Initialize(const double t0, const VecX &x0, const MatXd &P0);

    VecX GetMeasurementResidual() { return dy_; }
    MatX GetKalmanGain() { return K_; }
    MatX GetMeasurementNoiseCov() { return R_; }
    MatX GetMeasurementJacobian() { return H_; }
    int GetMeasurementSize() { return H_.rows(); }

    void SetOutlierThreshold(double outlier_threshold);
    int RemoveOutliers(int m);

    // Interface
    void SetDynamicsFunction(FilterDynamicsFunction f_dyn) { f_dyn_ = f_dyn; }
    void SetProcessNoiseFunction(FilterProcessNoiseFunction f_proc) { f_proc_ = f_proc; }
    void SetMeasurementFunction(FilterMeasurementFunction f_meas) { f_meas_ = f_meas; }

    void Predict(Real t_end);
    void Update(VecX z_true);

    VecX GetSate() { return x_; }
    VecX GetSatePrior() { return x_prior_; }
    VecX GetStatePost() { return x_post_; }

    MatXd GetCovariance() { return P_; }
    MatXd GetCovariancePrior() { return P_prior_; }
    MatXd GetCovariancePost() { return P_post_; }
  };
}  // namespace lupnt
