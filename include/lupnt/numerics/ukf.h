/**
 * @file ukf.h
 * @author Stanford NAV LAB
 * @brief  Unscented Kalman Filter
 * @version 0.1
 * @date 2024-12-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/numerics/filters.h"

namespace lupnt {

  class UKF : public IFilter {
  private:
    // Unscented transform parameters
    double alpha_ = 1e-3;
    double beta_ = 2.0;
    double kappa_ = 0.0;

    // Derived quantities
    int n_x_ = 0;        // dimension of the state
    int n_sigma_ = 0;    // number of sigma points = 2*n_x_ + 1
    double lambda_ = 0;  // alpha_^2*(n_x_+kappa_) - n_x_

    // Weights for mean and covariance
    VecXd w_m_;  // mean weights
    VecXd w_c_;  // covariance weights

    void InitializeUkfParams();

  protected:
    MatXd F_;  // State transition matrix
    MatXd H_;  // Measurement matrix
    MatXd Q_;  // Process noise cov
    MatXd R_;  // Measurement noise cov

    VecX dy_;       // Measurement residual
    VecX dx_;       // State update
    VecX z_true_;   // Observed measurement
    VecX z_prior_;  // Predicted measurement

    MatXd S_;  // Innovation cov
    MatXd K_;  // Kalman gain

  public:
    UKF() = default;

    /**
     * @brief Initialize the UKF filter epoch and state
     *
     * @param t0  Epoch time
     * @param x0  State vector
     * @param P0  Covariance matrix
     */
    void Initialize(const Real t0, const VecX& x0, const MatXd& P0) {
      t_ = t0;
      x_ = x0;
      P_ = P0;

      n_x_ = x0.size();
      x_prior_ = x0;
      P_prior_ = P0;
      x_post_ = x0;
      P_post_ = P0;

      InitializeUkfParams();
    };

    /**
     * @brief Compute the sigma points for the UKF
     *
     * @param state  State vector (n_x_)
     * @param cov  Covariance matrix (n_x_ x n_x_)
     * @return MatX  Sigma points (n_x_ x n_sigma_)
     */
    MatX ComputeSigmaPoints(const VecX& state, const MatXd& cov);

    /**
     * @brief Unscented transform
     *
     * @param sigma_points   Sigma points (n_x_ x n_sigma_)
     * @param state_mean     Mean of the state (n_x_)
     * @param cov_out        Covariance of the state (n_x_ x n_x_)
     * @param add_noise      Add process noise
     * @param Q              Process noise covariance (n_x_ x n_x_)
     */
    void UnscentedTransform(const MatX& sigma_points, VecX& state_mean, MatXd& cov_out,
                            bool add_noise, const MatXd& Q = MatXd::Zero(0, 0));

    void SetAlpha(double alpha) {
      alpha_ = alpha;
      InitializeUkfParams();
    }
    void SetBeta(double beta) {
      beta_ = beta;
      InitializeUkfParams();
    }
    void SetKappa(double kappa) {
      kappa_ = kappa;
      InitializeUkfParams();
    }

    void Predict(Real t_end) override;
    void Update(VecX z_true) override;

    // Interface
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