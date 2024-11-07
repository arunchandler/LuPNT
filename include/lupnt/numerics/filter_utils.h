/**
 * @file filter_utils.h
 * @author Stanford NAV Lab
 * @brief Utility functions for filters
 * @version 0.1
 * @date 2024-10-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/numerics/filters.h"
#include "lupnt/physics/clock.h"
#include "lupnt/agents/agent.h"

namespace lupnt {

  /**
   * @brief Construct the initial covariance matrix for the Orbit and Clock states
   *
   * @param pos_err initial position error [km]
   * @param vel_err initial velocity error [km/s]
   * @param clk_bias_err  initial clock bias error [s]
   * @param clk_drift_err  initial clock drift error [s/s]
   * @return MatXd  Initial covariance matrix
   */
  MatXd ConstructInitCovariancePVC(double pos_err, double vel_err, double clk_bias_err,
                                   double clk_drift_err);

  /**
   * @brief  Construct the process noise function for the RVC model
   *
   * @param cmodel  Clock model
   * @param state_size   State size
   * @param sigma_acc   Acceleration noise [km/s^2]
   * @return FilterProcessNoiseFunction
   */
  FilterProcessNoiseFunction ConstructProcessNoisePVC(ClockModel cmodel, int state_size,
                                                      double sigma_acc);

  /**
   * @brief Print the EKF progress header for the position, velocity, and clock states
   *
   */
  void PrintEKFProgressHeaderPVC();

  /**
   * @brief Compute the estimation error for the position, velocity, and clock states fpr the filter
   *
   * @param sat      Spacecraft
   * @param filter   Filter
   * @return VecXd   Estimation error vector (position, velocity, clock bias)
   */
  VecXd ComputeEstimationErrorPVC(const Ptr<Spacecraft> sat, IFilter* filter);

  /**
   * @brief Print the EKF progress for the position, velocity, and clock states
   *
   * @param t              Time
   * @param x_pos_err      Position error
   * @param x_vel_err      Velocity error
   * @param x_clk_bias_err Clock bias error
   */
  void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err);


}  // namespace lupnt
