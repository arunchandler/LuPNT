/**
 * @file gnssChannel.cpp
 * @author Stanford NAV LAB
 * @brief Gnss Channel
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/measurements/gnss_channel.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/measurements/gnss_receiver.h"
#include "lupnt/measurements/gnss_transmitter.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/spice_interface.h"

#define DEBUG_TRANSMISSIONS 0

namespace lupnt {
  /**
   * @brief Receiver receives all available Gnss signals
   *
   * @param rx
   * @param t
   * @return std::vector<GnssTransmission>
   */
  std::vector<GnssTransmission> GnssChannel::Receive(GnssReceiver &rx, Real t) {
    std::vector<GnssTransmission> received_transs;  // create an empty vector

    // Messages from other comms systems that can generate Gnss messages
    for (auto &tx : tx_devices) {
      // Solve light time delay
      Real tau = 0.0;  // light time delay
      CartesianOrbitState rv_rx_gcrf = rx.GetAgent()->GetCartesianGCRFStateAtEpoch(t);
      CartesianOrbitState rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);

      // Compute Light time delay
      Real tau_prev = 0.0;  // propagation time
      int max_iter = 100;
      for (int n_iter = 0; n_iter < max_iter; n_iter++) {
        rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
        Real rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm();
        tau = rho / C;
        if (abs(tau - tau_prev) < 1e-12)
          break;
        else {
          tau_prev = tau;
        }
      }

      // Transmission and reception times
      Real t_rx = t;
      Real t_tx = t - tau;

      // Convert to Moon Inertial frame
      auto rv_tx_mi = ConvertOrbitStateFrame(rv_tx_gcrf, t_tx, Frame::MOON_CI);
      auto rv_rx_mi = ConvertOrbitStateFrame(rv_rx_gcrf, t_rx, Frame::MOON_CI);

      // Occultation
      std::string tx_planet = "";
      std::map<std::string, bool> occult = Occultation::ComputeOccultationGnss(
          rv_tx_gcrf.r(), rv_tx_mi.r(), rv_rx_gcrf.r(), rv_rx_mi.r(), tx_planet, 10.0 * RAD);

      if (occult["earth"] || occult["moon"]) {
        // std::cout << "Earth or Moon occultation" << std::endl;
        continue;  // quit if occulted
      }

      Real Ar = rx.GetReceiverAntennaGain(t_rx, rv_tx_gcrf.r(), rv_rx_gcrf.r());

      // Generate transmission
      GnssTransmission trans = tx->GenerateTransmission(t_tx);
      Real d = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm();

      // Link budget
      for (size_t freq_idx = 0; freq_idx < tx->freq_list.size(); freq_idx++) {
        std::string freq_name = tx->freq_list[freq_idx];
        Real freq = tx->freq_map[freq_name];
        Real Ad = 20.0 * log10((C / freq) / (4.0 * PI * d));
        Real scalars = tx->P_tx + rx.rx_param_.Ae + rx.rx_param_.As
                       - (10.0 * log10(rx.rx_param_.Tsys)) + 228.6 + rx.rx_param_.L;

        // Transmitter and Receiver Antenna gain
        Real At
            = tx->GetTransmitterAntennaGainFreq(t_tx, rv_tx_gcrf.r(), rv_rx_gcrf.r(), freq_name);
        trans.CN0 = At + Ar + Ad + scalars;

        if (std::isnan(At.val()) || occult["earth"] || occult["moon"]
            || trans.CN0 < rx.rx_param_.CN0threshold) {
          // not visible
          continue;
        } else {
          trans.AP = tx->P_tx + At + Ad + rx.rx_param_.Ae;
          trans.RP = trans.AP + Ar + rx.rx_param_.As;
          trans.vis_antenna = true;

          // TX
          trans.t_tx = t_tx;
          trans.freq = freq;
          trans.freq_label = freq_name;
          trans.chip_rate = tx->rc_map[freq_name];
          trans.dt_tx = 0.0;      // Todo: Get this from ephemeris
          trans.dt_tx_dot = 0.0;  // Todo: Get this from ephemeris
          trans.r_tx = rv_tx_gcrf.r().cast<double>();
          trans.v_tx = rv_tx_gcrf.v().cast<double>();

          // Channel
          trans.I_rx = 0.0;
          trans.T_rx = 0.0;
          trans.vis_atmos = 1 - occult["atmos"];
          trans.vis_ionos = 1 - occult["ionos"];
          trans.vis_earth = 1 - occult["earth"];
          trans.vis_moon = 1 - occult["moon"];

          trans.ID_tx = tx->GetPRN();

          // RX
          trans.t_rx = t_rx;
          trans.dt_rx = rx.GetAgent()->GetClockState().GetValue(0).val();
          trans.dt_rx_dot = rx.GetAgent()->GetClockState().GetValue(1).val();
          trans.r_rx = rv_rx_gcrf.r().cast<double>();
          trans.v_rx = rv_rx_gcrf.v().cast<double>();

          // receiver chip param
          trans.gnssr_param = rx.gnssr_param_;

          // std::cout << "prn: " << tx->GetPRN() << " freq: " << freq_name << " At:" << At << "
          // Ar:" << Ar << "  C/N0:" << trans.CN0 << std::endl;

          received_transs.push_back(trans);
        }
      }
    }
    return received_transs;
  }  // Receive

}  // namespace lupnt
