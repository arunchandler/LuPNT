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

#include "lupnt/agents/spacecraft.h"
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
  MatXd InitialCovariancePosVelClock(double pos_err, double vel_err, double clk_bias_err,
                                     double clk_drift_err);

  /**
   * @brief Construct the process noise function for the linear RV model
   *
   * @param sigma_acc  Acceleration noise [km/s^2]
   * @return FilterProcessNoiseFunction
   */
  FilterProcessNoiseFunction ProcessNoiseLinearPosVel(double sigma_acc);

  /**
   * @brief Construct the process noise function for the clock model
   *
   * @param cmodel  Clock model
   * @param clock_state_size   Clock state size
   * @return FilterProcessNoiseFunction
   */
  FilterProcessNoiseFunction ProcessNoiseClock(ClockModel cmodel, int clock_state_size);

  /**
   * @brief Construct a new Filter Process Noise Function P V C object for the EKF
   *
   * @param cmodel
   * @param state_size
   * @param sigma_acc
   * @param n_sat
   */
  FilterProcessNoiseFunction ProcessNoisePosVelClock(ClockModel cmodel, int state_size,
                                                     double sigma_acc, int n_sat = 1);

  /**
   * @brief Print the EKF progress header for the position, velocity, and clock states
   *
   */
  void PrintEKFProgressHeaderPVC();

  /**
   * @brief Print the EKF progress header for the position, velocity, and clock states for N
   * satellites
   *
   * @param n_sat  Number of satellites
   */
  void PrintEKFProgressHeaderPVC(int n_sat);

  /**
   * @brief Construct the true state vector from the spacecraft vector
   *
   * @param sats <std::vector<Ptr<Spacecraft>>  Spacecraft vector
   * @return VecX  True state vector (8xN)
   */
  VecX ConstructTrueStateVecFromSats(const std::vector<Ptr<Spacecraft>>& sats);

  /**
   * @brief Compute the estimation error for the position, velocity, and clock states fpr the filter
   *
   * @param sat      Spacecraft
   * @param filter   Filter
   * @param start_idx  Start index for the filter state vector
   * @return VecXd   Estimation error vector (position, velocity, clock bias)
   */
  VecXd ComputeEstimationErrorPVC(const Ptr<Spacecraft> sat, IFilter* filter, int start_idx);

  /**
   * @brief Compute the estimation error for the position, velocity, and clock states for the filter
   *
   * @param sats    Spacecraft vector
   * @param filter  Filter
   * @return VecXd  Estimation error vector (position, velocity, clock bias)
   */
  VecXd ComputeEstimationErrorPVC(const std::vector<Ptr<Spacecraft>>& sats, IFilter* filter);

  /**
   * @brief Print the EKF progress for the position, velocity, and clock states
   *
   * @param t              Time
   * @param x_pos_err      Position error
   * @param x_vel_err      Velocity error
   * @param x_clk_bias_err Clock bias error
   */
  void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err);

  /**
   * @brief Print the EKF progress for the position, velocity, and clock states for N satellites
   *
   * @param t              Time
   * @param est_err        Estimation error vector (position, velocity, clock bias)
   * @param n_sat          Number of satellites
   *
   */
  void PrintEKFProgressPVC(double t, const VecXd& est_err, int n_sat);

  /**
   * @brief Print the State Estimation Error Statistics
   *
   * @param num_meas  Number of GPS measurements
   * @param error_mat Error matrix (position, velocity, clock bias, clock drift)
   * @param ratio     Ratio
   * @param n_sat     Number of satellites
   *
   */
  void PrintEstimationStatistics(const VecXd& num_meas, const MatXd& error_mat, double ratio,
                                 int n_sat);

}  // namespace lupnt
