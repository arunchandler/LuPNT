/**
 * @file RadioMeasurement.cpp
 * @author Stanford NAV LAB
 * @brief Class for Radionavigation measurements
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "lupnt/agents/agent.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  /**
   * @brief Compute the one-way range between two points
   *
   * @param r_tx   Transmitter position
   * @param r_rx   Receiver position
   * @param offset  Measurement offset
   * @return Real   One-way range
   */
  Real ComputeOneWayRange(const VecX& r_tx, const VecX& r_rx, Real offset);

  /**
   * @brief Compute the pseudorange between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param dt_tx  Transmitter clock offset
   * @param dt_rx  Receiver clock offset
   * @param offset  Measurement offset
   * @return Real  Pseudorange
   */
  Real ComputePseudorange(const VecX& r_tx, const VecX& r_rx, Real dt_tx, Real dt_rx, Real offset);

  /**
   * @brief Compute the pseudorangerate between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param v_tx  Transmitter velocity
   * @param v_rx  Receiver velocity
   * @param dt_tx_dot  Transmitter clock offset rate
   * @param dt_rx_dot  Receiver clock offset rate
   * @param offset  Measurement offset
   * @return Real   Pseudorangerate
   */
  Real ComputePseudorangerate(const VecX& r_tx, const VecX& r_rx, const VecX& v_tx,
                              const VecX& v_rx, Real dt_tx_dot, Real dt_rx_dot, Real offset);

  /**
   * @brief Compute the Doppler shift between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param v_tx  Transmitter velocity
   * @param v_rx  Receiver velocity
   * @param dt_tx_dot  Transmitter clock offset rate
   * @param dt_rx_dot  Receiver clock offset rate
   * @param offset  Measurement offset
   * @return Real  Doppler shift [Hz]
   */
  Real ComputeDopplerShift(const VecX& r_tx, const VecX& r_rx, const VecX& v_tx, const VecX& v_rx,
                           Real dt_tx_dot, Real dt_rx_dot, Real f, Real offset);

  /**
   * @brief Compute the one-way range between two points considering light time
   * delay
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param epoch_rx_recorded   Reception epoch (TAI) t_R
   * @param epoch_ref  Reference epoch
   * @param rv_tx  Transmitter position at epoch_ref (w.r.t to
   * central body)
   * @param rv_rx   Receiver position at epoch_ref  (w.r.t to
   * central body)
   * @param dt_tx      Transmitter clock offset and bias at epoch_ref
   * @param dt_rx      Receiver clock offset and bias at epoch_ref
   * @param agent_tx  Transmitter agent
   * @param agent_rx  Receiver agent
   * @param additional_delay  Additional delay for the receiver time recording
   * @param use_pure_range  Flag to indicate if the pure range is used (no clock bias)
   * @return Real      One-way pseudorange at t_R (clock offset error included)
   */
  Real ComputeOneWayRangeLTR(Real epoch_rx_recorded, Real epoch_ref, const Vec6& rv_tx,
                             const Vec6& rv_rx, const Vec2& dt_tx, const Vec2& dt_rx,
                             Ptr<Agent> agent_tx, Ptr<Agent> agent_rx, Real addional_delay,
                             bool use_pure_range);

  /**
   * @brief Compute the one-way range between two points considering light time
   * delay
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param epoch_rx_recorded   Reception epoch (TAI) t_R
   * @param epoch_ref  Reference epoch
   * @param rv_tx  Transmitter position at epoch_ref (w.r.t to
   * central body)
   * @param rv_rx   Receiver position at epoch_ref  (w.r.t to
   * central body)
   * @param dt_tx      Transmitter clock offset and bias at epoch_ref
   * @param dt_rx      Receiver clock offset and bias at epoch_ref
   * @param agent_tx  Transmitter agent
   * @param agent_rx  Receiver agent
   * @param hardware_delay_target  Hardware delay at the target
   * @param additional_delay Additional delay for the receiver time
   * @param use_pure_range  Flag to indicate if the pure range is used (no clock bias)
   * @return Real      One-way pseudorange at t_R (clock offset error included)
   */
  Real ComputeTwoWayRangeLTR(Real epoch_rx, Real epoch_ref, const Vec6& rv_target,
                             const Vec6& rv_rx, const Vec2& clk_target, const Vec2& clk_receiver,
                             Ptr<Agent> agent_target, Ptr<Agent> agent_rx,
                             Real hardware_delay_target, Real additional_delay = 0.0);

  /**
   * @brief Compute the one-way range rate between two points considering light
   * time delay
   *
   * @param epoch_rx Reception epoch (TAI) t_R
   * @param rv_tx_tr  Transmitter position at reception time t_R (w.r.t to
   * @param rv_rx_tr  Receiver position at reception time t_R (w.r.t to
   * @param dt_dot_tx  Transmitter clock offset rate
   * @param dt_dot_rx  Receiver clock offset rate
   * @param target_center_body  Central body of the target
   * @param rx_center_body  Central body of the receiver
   * @param is_bodyfixed_target  Flag to indicate if the target state is
   * body-fixed
   * @param is_bodyfixed_rx  Flag to indicate if the receiver state is
   * body-fixed
   * @param hardware_delay  Hardware delay [s]
   * @param T_I  Integration time [s]
   * @return Real
   */
  Real ComputeOneWayRangeRateLTR(Real epoch_rx, Real epoch_ref, const Vec6& rv_tx_tr,
                                 const Vec6& rv_rx_tr, const Vec2& clk_tx, const Vec2& clk_rx,
                                 Ptr<Agent> agent_tx, Ptr<Agent> agent_rx, Real hardware_delay,
                                 Real T_I);

  /**
   * @brief Compute the two-way range rate between two points considering light
   *
   * @param epoch_rx        Reception epoch (TAI) t_R
   * @param epoch_ref       Reference epoch
   * @param rv_target_tr    Transmitter position at reception time t_R (w.r.t to
   * @param rv_rx_tr        Receiver position at reception time t_R (w.r.t to
   * @param clk_target      Transmitter clock state at the ref epoch
   * @param clk_rx          Receiver clock state at the ref epoch
   * @param agent_target    Transmitter agent
   * @param agent_receiver  Receiver agent
   * @param hardware_delay  Hardware delay at the target relay
   * @param T_I             Integration time [s]
   * @return Real
   */
  Real ComputeTwoWayRangeRateLTR(Real epoch_rx, Real epoch_ref, const Vec6& rv_target_tr,
                                 const Vec6& rv_rx_tr, const Vec2& clk_target, const Vec2& clk_rx,
                                 Ptr<Agent> agent_target, Ptr<Agent> agent_receiver,
                                 Real hardware_delay, Real T_I);

  /**
   * @brief Compute the PN regenerative range error for chip tracking loop
   * (CTL): suited for onboard processing
   *  Reference: "Pseudo-Noise (PN) Ranging Systems", Greenbook 2014
   *
   * @param PRC_N0   Carrier-to-noise ratio for the range clock (linear)
   * @param B_L_CTL      One-sided Chip tracking Loop noise bandwidth
   * (usually around 1Hz, 0.5Hz)
   * @param T_c       Chip period [s] = = 1/ (2 f_RC)
   * @return Real range error [m]
   */
  Real ComputePnRangeErrorCTL(Real PRC_N0, Real B_L_CTL, Real T_c,
                              Modulation modulation_type = Modulation::BPSK);

  /**
   *
   * @brief Compute the PN regenerative range error for open loop (OL) tracking
   * suited for ground stations
   * Reference: "Pseudo-Noise (PN) Ranging Systems", Greenbook 2014
   *
   * @param PRC_N0 Carrier-to-noise ratio for the range clock (linear)
   * @param T_I Integration time
   * @param T_c Chip period [s] = 1/ (2 f_RC)
   * @return Real   range error [m]
   */
  Real ComputePnRangeErrorOL(Real PRC_N0, Real T_I, Real T_c,
                             Modulation modulation_type = Modulation::BPSK);

  /**
   * @brief Compute range rate error
   * Reference: https://deepspace.jpl.nasa.gov/dsndocs/810-005/202/202E.pdf
   *
   * @param B_L_carrier
   * @param f_C   downlink carrier frequency  [Hz]
   * @param T_s   period of the binary symbol [s]
   * @param T_I   integration time [s]
   * @param sigma_y_1s  Allan deviation of the receiver clock at 1s
   * @param PT_N0  [Hz]
   * @return Real  range rate error [m/s]
   */
  Real ComputeRangeRateErrorOneWay(Real B_L_carrier, Real f_C, Real T_s, Real T_I, Real PT_N0,
                                   Real sigma_y_1s, Modulation modulation_type = Modulation::BPSK,
                                   Real m_R = 0.0);

  /**
   * @brief Compute range rate error for two-way ranging
   * Reference: https://deepspace.jpl.nasa.gov/dsndocs/810-005/202/202E.pdf
   *
   * @param B_L_carrier  carrier loop noise bandwidth [Hz]
   * @param f_C   downlink carrier frequency [Hz]
   * @param T_s   period of the binary symbol [s]
   * @param T_I   integration time [s]
   * @param sigma_y_1s  Allan deviation of the receiver clock at 1s
   * @param G     transponder turnaround ratio
   * @param PT_N0  downlink total signal power to noise spectral density ratio
   * [Hz]
   * @return Real  range rate error [m/s]
   */
  Real ComputeRangeRateErrorTwoWay(Real B_L_carrier, Real f_C, Real T_s, Real T_I, Real PT_N0,
                                   Real sigma_y_1s, Real G,
                                   Modulation modulation_type = Modulation::BPSK, Real m_R = 0.0);

}  // namespace lupnt
