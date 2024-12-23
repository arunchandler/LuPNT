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
  public:
    virtual ~IFilter() = default;

    virtual void SetDynamicsFunction(FilterDynamicsFunction f_dyn) = 0;
    virtual void SetProcessNoiseFunction(FilterProcessNoiseFunction f_proc) = 0;
    virtual void SetMeasurementFunction(FilterMeasurementFunction f_meas) = 0;

    virtual void Predict(Real t_end) = 0;
    virtual void Update(VecX z_obs) = 0;

    virtual VecX GetSate() = 0;
    virtual VecX GetSatePrior() = 0;
    virtual VecX GetStatePost() = 0;
    virtual MatXd GetCovariance() = 0;
    virtual MatXd GetCovariancePrior() = 0;
    virtual MatXd GetCovariancePost() = 0;
  };

}  // namespace lupnt
