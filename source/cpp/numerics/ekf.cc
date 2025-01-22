/**
 * @file filters.cpp
 * @author Stanford NAVLAB
 * @brief Implementation of Filters
 * @version 0.1
 * @date 2023-09-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/numerics/ekf.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  /*****************************************************
   *   Extended Kalman Filter
   *****************************************************/

  void EKF::Initialize(const Real t0, const VecX &x0, const MatXd &P0) {
    t_ = t0;
    x_ = x0;
    P_ = P0;
    P_prior_ = P0;
    x_prior_ = x0;
    P_post_ = P0;
    x_post_ = x0;

    // Set Matrix
    F_.resize(x_.size(), x_.size());
  }

  void EKF::Predict(Real t_end) {
    if ((adaptive_process_noise_) && ((Q_set_) && (S_set_))) {
      MatXd Q_tilde = dx_.cast<double>() * dx_.transpose().cast<double>();
      Q_ = alpha_Q_ * Q_ + (1 - alpha_Q_) * Q_tilde;
    } else {
      Q_ = f_proc_(x_, t_, t_end);
      Q_set_ = true;
    }

    x_ = f_dyn_(x_, t_, t_end, &F_);

    P_ = F_ * P_ * F_.transpose() + Q_;
    x_prior_ = x_;
    P_prior_ = P_;
    t_ = t_end;
  }

  void EKF::SetOutlierThreshold(double outlier_threshold) {
    if (outlier_threshold < 0) throw std::invalid_argument("Outlier threshold must be positive");
    outlier_threshold_ = outlier_threshold;
  }

  /**
   * @brief Update step
   *
   * @param z_true observed measurement
   */
  void EKF::Update(VecX z_true_in) {
    z_true_ = z_true_in;
    x_post_ = x_;
    P_post_ = P_;

    int n = x_.size();
    int m = z_true_.size();
    if (m == 0) return;  // no measurement, nothing to update

    // Prior measurement
    z_prior_ = f_meas_(x_, &H_, &R_);

    S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information
    S_set_ = true;

    dy_ = z_true_ - z_prior_;

    // Remove outliers
    m = RemoveOutliers(m);
    if (m == 0) return;  // all measurements are outliers

    // Update step
    MatXd S_inv = S_.inverse();
    // MatXd S_inv = S_.completeOrthogonalDecomposition().pseudoInverse();
    K_ = P_ * H_.transpose() * S_inv;  // Kalman gain

    // K_ = P_ * H_.transpose() * PseudoInverse(S_);

    dx_ = K_ * dy_;
    x_ = x_ + dx_;
    MatXd I = MatXd::Identity(n, n);
    MatXd G(n, n);
    G = I - K_ * H_;

    P_ = G * P_ * G.transpose() + K_ * R_ * K_.transpose();  // Joseph form

    // covariance inflation
    double lambda = 0.0;
    P_ = (1 + lambda) * P_;
    // P_ = G * P_;

    x_post_ = x_;
    P_post_ = P_;
  }

  /**
   * @brief Remove outliers from the measurement
   *
   * @param m   number of measurements
   * @param debug   debug flag
   * @return int   number of measurements after removing outliers
   */
  int EKF::RemoveOutliers(int m_orig) {
    std::vector<int> is_outlier(m_orig);
    VecXd ratio(m_orig);
    int n_valid = 0;
    int n = x_.size();

    for (int i = 0; i < m_orig; i++) {
      ratio(i) = abs(dy_(i).val() / sqrt(S_(i, i)));
      if (ratio(i) > outlier_threshold_) {
        is_outlier[i] = 1;
      } else {
        is_outlier[i] = 0;
        n_valid++;
      }
    }

    // Remove outliers
    int m = n_valid;
    if (m == m_orig) {
      return m;  // all measurement valid, nothing to change
    }

    MatXd H_new(m, n);
    MatXd R_new(m, m);
    MatXd S_new(m, m);
    VecX dy_new(m);
    int j = 0;
    int l = 0;

    for (int i = 0; i < m_orig; i++) {
      if (is_outlier[i] == 0) {
        H_new.row(j) = H_.row(i);
        dy_new(j) = dy_(i);

        l = 0;
        for (int k = 0; k < m_orig; k++) {
          if (is_outlier[k] == 0) {
            R_new(j, l) = R_(i, k);
            S_new(j, l) = S_(i, k);
            l++;
          }
        }

        j++;
      }
    }

    // set new value
    H_ = H_new;
    R_ = R_new;
    S_ = S_new;
    dy_ = dy_new;

    return m;
  }

}  // namespace lupnt
