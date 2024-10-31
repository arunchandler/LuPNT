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
  MatXd ConstructInitCovarianceRVC(double pos_err, double vel_err, double clk_bias_err,
                                   double clk_drift_err);

  /**
   * @brief  Construct the process noise function for the RVC model
   *
   * @param cmodel  Clock model
   * @param state_size   State size
   * @param sigma_acc   Acceleration noise [km/s^2]
   * @return FilterProcessNoiseFunction
   */
  FilterProcessNoiseFunction ConstructProcessNoiseRVC(ClockModel cmodel, int state_size,
                                                      double sigma_acc);

}  // namespace lupnt
