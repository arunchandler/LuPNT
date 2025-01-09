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
  protected:
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

    bool adaptive_process_noise_ = false;
    bool Q_set_ = false;
    bool S_set_ = false;

    // adaptive process noise parameters
    double alpha_Q_ = 0.95;

  public:
    double outlier_threshold_ = 3.0;

    void Initialize(const double t0, const VecX &x0, const MatXd &P0);

    void SetOutlierThreshold(double outlier_threshold);
    int RemoveOutliers(int m);

    void SetAdaptiveProcessNoise(bool adaptive_process_noise) {
      adaptive_process_noise_ = adaptive_process_noise;
      Q_set_ = false;
      S_set_ = false;
    }

    void SetAdaptiveProcessNoiseCoeff(double alpha_Q) { alpha_Q_ = alpha_Q; }

    // Interface
    void Predict(Real t_end) override;
    void Update(VecX z_true) override;

    VecXd GetMeasurementResidual() { return dy_.cast<double>(); }
    MatXd GetKalmanGain() { return K_; }
    MatXd GetMeasurementNoiseCov() { return R_; }
    MatXd GetMeasurementJacobian() { return H_; }
    int GetMeasurementSize() { return H_.rows(); }
    MatXd GetProcessNoise() { return Q_; }
    MatXd GetStateJacobian() { return F_; }
    MatXd GetInnovationCov() { return S_; }
    MatXd GetMeasurementCov() { return R_; }
    VecXd GetStateCorrection() { return dx_.cast<double>(); }
    VecXd GetTrueMeasurement() { return z_true_.cast<double>(); }
    VecXd GetPredictedMeasurement() { return z_prior_.cast<double>(); }
  };
}  // namespace lupnt
