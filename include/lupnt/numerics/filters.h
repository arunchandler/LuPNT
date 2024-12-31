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

namespace lupnt {

  // Dynamics and Measurement Function

  /**
   * @brief Dynamics function for the filter
   *
   * @param x State
   * @param t_curr Current time
   * @param t_end End time
   * @param Phi STM of the dynamics
   */
  typedef std::function<VecX(const VecX, Real t_curr, Real t_end, MatXd *F)> FilterDynamicsFunction;

  /**
   * @brief Process noise function for the filter
   *
   * @param x State
   * @param t_curr Current time
   * @param t_end End time
   * @return VecXd Process noise covariance
   *
   */
  typedef std::function<MatXd(const VecX, Real t_curr, Real t_end)> FilterProcessNoiseFunction;

  /**
   * @brief Measurement function for the filter
   *
   * @param x State
   * @param H Measurement matrix
   * @param R Measurement noise covariance
   *
   */
  typedef std::function<VecX(const VecX x, MatXd *H, MatXd *R)> FilterMeasurementFunction;

  class IFilter {

  protected:
    FilterDynamicsFunction f_dyn_;
    FilterProcessNoiseFunction f_proc_;
    FilterMeasurementFunction f_meas_;

    Real t_;
    VecX x_;
    MatXd P_;
    VecX x_prior_;
    VecX x_post_;
    MatXd P_prior_;
    MatXd P_post_;

  public:
    virtual ~IFilter() = default;

    void SetDynamicsFunction(FilterDynamicsFunction f_dyn) { f_dyn_ = f_dyn; }
    void SetProcessNoiseFunction(FilterProcessNoiseFunction f_proc) { f_proc_ = f_proc; }
    void SetMeasurementFunction(FilterMeasurementFunction f_meas) { f_meas_ = f_meas; }

    virtual void Predict(Real t_end) = 0;
    virtual void Update(VecX z_obs) = 0;

    VecX GetState() { return x_; }
    VecX GetStatePrior() { return x_prior_; }
    VecX GetStatePost() { return x_post_; }

    MatXd GetCovariance() { return P_; }
    MatXd GetCovariancePrior() { return P_prior_; }
    MatXd GetCovariancePost() { return P_post_; }
  };

}  // namespace lupnt
