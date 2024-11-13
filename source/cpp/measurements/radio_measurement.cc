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

#include "lupnt/measurements/radio_measurement.h"

#include "lupnt/measurements/comm_utils.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {
  Real ComputeOneWayRange(VecX r_tx, VecX r_rx, Real offset) {
    Real rho_rx = (r_tx - r_rx).norm();
    return rho_rx + offset;
  };

  Real ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx, Real offset) {
    // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
    Real rho_rx = (r_tx - r_rx).norm();
    Real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
    return P_rx;
  };

  Real ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                              Real dt_rx_dot, Real offset) {
    VecX e_rx = (r_tx - r_rx).normalized();
    Real prr = e_rx.dot(v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
    return prr;
  };

  Real ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                           Real dt_rx_dot, Real f, Real offset) {
    Real f_D
        = -f / C * ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot, dt_rx_dot, offset);
    return f_D;
  };

  Real ComputeOneWayRangeLTR(Real epoch_rx_local, Real epoch_ref, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                             const Ptr<Agent> agent_tx, const Ptr<Agent> agent_rx, Real additional_delay, bool use_pure_range=false) {

    // solve for tau_d (downlink time)
    int max_iter = 5;
    Real tau_d = 0.0;   // light time
    Real tau_d_prev = 0.0;

    Vec3 r0_p, rid_p, rho_ad;
    Real r0_p_norm, rid_p_norm, rho_ad_norm;
    Real rx_clk_offset = 0.0;
    Real rx_epoch = epoch_rx_local;
    Real tx_epoch = epoch_rx_local;

    // Solve for Downlink
    for (int i = 0; i < max_iter; i++) {
      if (epoch_ref == epoch_rx_local) {
        rx_clk_offset = clk_rx(0);
      }
      else {
        rx_clk_offset = agent_rx->PropagateClockState(epoch_ref, clk_rx, epoch_rx_local)(0);
      }
      rx_epoch = epoch_rx_local - rx_clk_offset - additional_delay;
      tx_epoch = rx_epoch - tau_d;

      // Propagate the agent from the current epoch to the downlink/uplink  epoch
      r0_p = agent_rx->PropagateRvState(epoch_ref, rv_rx, rx_epoch, Frame::GCRF).segment(0, 3);
      rid_p = agent_tx->PropagateRvState(epoch_ref, rv_tx, tx_epoch, Frame::GCRF).segment(0, 3);
      rho_ad = rid_p - r0_p;

      // norms
      r0_p_norm = r0_p.norm();
      rid_p_norm = rid_p.norm();
      rho_ad_norm = rho_ad.norm();

      // Compensate for relativistic effects (Shapiro time delay)
      // double shapiro = 2 * mu_rx / C *
      //                  log((r0_p_norm + rid_p_norm + rho_ad_norm) /
      //                      (r0_p_norm + rid_p_norm - rho_ad_norm));

      tau_d = rho_ad.norm() / C;  // + shapiro;

      if (fabs(tau_d.val() - tau_d_prev.val()) < 1e-13) {  // goes under pico-second
        break;
      } else {
        tau_d_prev = tau_d;
      }
    }

    Real tx_clk_offset = agent_tx->PropagateClockState(epoch_ref, clk_tx, tx_epoch)(0);
    Real epoch_tx_local = tx_epoch + tx_clk_offset;

    Real rho_d = 0.0;
    if (use_pure_range) {
      rho_d = rho_ad_norm;
    }
    else {
      rho_d = C * (epoch_rx_local - epoch_tx_local);
    }

    return rho_d;
  };

  Real ComputeTwoWayRangeLTR(Real epoch_rx_local, Real epoch_ref, Vec6 rv_target, Vec6 rv_rx,
                             Vec2 clk_target, Vec2 clk_receiver,
                             Ptr<Agent> agent_target, Ptr<Agent> agent_receiver, Real hardware_delay_target,
                             Real additional_delay) {

    // solve for tau_d (downlink time, target->rx)
    Real rho_d = ComputeOneWayRangeLTR(epoch_rx_local, epoch_ref, rv_target, rv_rx, clk_target, clk_receiver,
                                       agent_target, agent_receiver, additional_delay, true);
    Real tau_d = rho_d / C;

    // solve for tau_u (uplink time, rx->target)
    Real delay_uplink = additional_delay + tau_d + hardware_delay_target;  // total delay for uplink w.r.t epoch_rx
    Real rho_u = ComputeOneWayRangeLTR(epoch_rx_local, epoch_ref, rv_rx, rv_target, clk_receiver, clk_target,
                                       agent_receiver, agent_target, delay_uplink, true);
    Real tau_u = rho_u / C;

    Real rho_ud = C / 2 * (tau_u + tau_d);

    return rho_ud;
  };

  Real ComputeOneWayRangeRateLTR(Real epoch_rx, Real epoch_ref, Vec6 rv_tx_tr, Vec6 rv_rx_tr, Vec2 clk_tx, Vec2 clk_rx, 
                                 Ptr<Agent> agent_tx, Ptr<Agent> agent_rx, Real additional_delay, double T_I) {
    Real rho_d = ComputeOneWayRangeLTR(epoch_rx, epoch_ref, rv_tx_tr, rv_rx_tr, clk_tx, clk_rx,
                                       agent_tx, agent_rx, additional_delay, false);
    Real rho_d_past = ComputeOneWayRangeLTR(epoch_rx, epoch_ref, rv_tx_tr, rv_rx_tr, clk_tx, clk_rx,
                                            agent_tx, agent_rx, additional_delay + T_I, false);

    Real rho_dot = (rho_d - rho_d_past) / T_I;

    return rho_dot;
  }

  Real ComputeTwoWayRangeRateLTR(Real epoch_rx, Real epoch_ref, Vec6 rv_target_tr, Vec6 rv_rx_tr, Vec2 clk_target, Vec2 clk_receiver,
                                 Ptr<Agent> agent_target, Ptr<Agent> agent_receiver, Real hardware_delay, double T_I) {
    Real rho_ud = ComputeTwoWayRangeLTR(epoch_rx, epoch_ref, rv_target_tr, rv_rx_tr, clk_target, clk_receiver, agent_target, agent_receiver,
                                        hardware_delay, 0);
    Real rho_ud_past = ComputeTwoWayRangeLTR(epoch_rx, epoch_ref, rv_target_tr, rv_rx_tr, clk_target, clk_receiver, agent_target, agent_receiver,
                                             hardware_delay, T_I);

    Real rho_dot = (rho_ud - rho_ud_past) / T_I;

    return rho_dot;
  }

  double ComputePnRangeErrorCTL(double PRC_N0, double B_L, double Tc, Modulation modulation_type) {
    (void)modulation_type;
    double sigma = 0.0;
    double f_RC = 1 / (2 * Tc);

    // Thermal noise
    sigma = 1 / sqrt(2) * C / (8 * f_RC) * sqrt(B_L / PRC_N0);

    // Error degrade for GMSK + PN

    return sigma;
  }

  double ComputePnRangeErrorOL(double PRC_N0, double TI, double Tc, Modulation modulation_type) {
    (void)modulation_type;
    double sigma = 0.0;
    double f_RC = 1 / (2 * Tc);

    // Thermal noise
    sigma = 1 / sqrt(32 * PI * PI) * (C / f_RC) * sqrt(1 / PRC_N0 / TI);

    // Error degrade for GMSK + PN

    return sigma;
  }

  double ComputeRangeRateErrorOneWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s, Modulation modulation_type,
                                     double m_R) {
    // Thermal noise
    double rho_L = ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, modulation_type, m_R);
    double sigma_vn = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I);

    // phase noise contribution
    double sigma_y_T = sigma_y_1s / sqrt(T_I);
    double sigma_vf = C * sigma_y_T;

    // phase scintillation
    double sigma_vs = 0.0;  // Asssume 0

    // Total Doppler Error
    double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));

    return sigma_v;
  }

  double ComputeRangeRateErrorTwoWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s, double G,
                                     Modulation modulation_type, double m_R) {
    // Thermal noise
    double rho_L = ComputeCarrierLoopSNR(PT_N0, B_L_carrier, T_s, modulation_type, m_R);
    double sigma_vnu = sqrt(1 / 2) * (C / (2 * PI * f_C * T_I)) * G / sqrt(rho_L);
    double sigma_vnd = sqrt(2 / rho_L) * C / (2 * PI * f_C * T_I) / sqrt(rho_L);

    double sigma_vn = sqrt(pow(sigma_vnu, 2) + pow(sigma_vnd, 2));

    // phase noise contribution
    double sigma_y_T = sigma_y_1s / sqrt(T_I);
    double sigma_vf = C * sigma_y_T / sqrt(2);

    // phase scintillation
    double sigma_vs = 0.0;  // Asssume 0

    // Total Doppler Error
    double sigma_v = sqrt(pow(sigma_vn, 2) + pow(sigma_vf, 2) + pow(sigma_vs, 2));

    // Error degrade for GMSK + PN

    return sigma_v;
  }

}  // namespace lupnt
