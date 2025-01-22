/**
 * @file SpaceChannel.cpp
 * @author Stanford NAV LAB
 * @brief Base spacechannel class (Under devlopment )
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/space_channel.h"

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/transmission.h"
#include "lupnt/physics/occultation.h"

namespace lupnt {

  ITransmission SpaceChannel::ComputeLinkBudget(Transmitter *tx, Receiver *rx, Real t,
                                                std::string time_fixed, bool compute_cn0) {
    ITransmission trans;  // create an empty vector

    // Transmitter and receiver positions and velocities
    Real tau = 0.0;  // light time delay
    CartesianOrbitState rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    CartesianOrbitState rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    Real t_rx = t;
    Real t_tx = t;

    // Solve light time delay
    if (time_fixed == "rx") {
      tau = SolveLightTimeDelayRx(tx, rx, t);
      rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
      t_rx = t;
      t_tx = t - tau;
    } else if (time_fixed == "tx") {
      tau = SolveLightTimeDelayTx(tx, rx, t);
      rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
      t_rx = t + tau;
      t_tx = t;
    } else {
      // generate error
      std::cerr << "Error: Invalid time_fixed parameter" << std::endl;
    }

    // assign values to trans
    trans.t_tx = t_tx.val();
    trans.t_rx = t_rx.val();
    trans.r_tx = rv_tx_gcrf.r().cast<double>();
    trans.v_tx = rv_tx_gcrf.v().cast<double>();
    trans.r_rx = rv_rx_gcrf.r().cast<double>();
    trans.v_rx = rv_rx_gcrf.v().cast<double>();

    // Commpute Occultations
    bool vis_all = true;
    std::map<std::string, bool> vis_occult;

    if (occult_bodies_.size() > 0) {
      Real epoch = (t_tx + t_rx) / 2.0;
      vis_occult = Occultation::ComputeOccultation(
          epoch, rv_tx_gcrf.r(), rv_rx_gcrf.r(), Frame::GCRF, Frame::GCRF, occult_bodies_,
          occult_alt_, elev_masks_, use_elev_mask_tx_, use_elev_mask_rx_);
      vis_all = vis_occult["all"];
    }

    // Link Budget
    if (compute_cn0) {
      Real At = tx->GetTransmitterAntennaGain(t_tx.val(), trans.r_tx, trans.r_rx);
      Real Ar = rx->GetReceiverAntennaGain(t_rx.val(), trans.r_tx, trans.r_rx);

      Real dist = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm();
      Real lambda = C / tx->freq_tx;
      Real fsl_loss_dB = ComputeFreeSpaceLossdB(dist, lambda);

      Real EIRP_dB = tx->P_tx + At;
      Real G_T_rx_dB = Ar - 10.0 * log10(rx->rx_param_.Tsys);
      Real loss = rx->rx_param_.Ae + rx->rx_param_.As + rx->rx_param_.L;  // sum of lossess (minus)
      Real CN0 = EIRP_dB - fsl_loss_dB + 228.6 + G_T_rx_dB + loss;
      // double CN = CN0 - 10.0 * log10(tx->bandwidth);
      // double RP = CN0 - 228.6 + 10.0 * log10(rx->rx_param_.Tsys);  // Received Power
      // double RP_N0 = RP - 10.0 * log10(rx->rx_param_.Tsys);        // Received Power Noise

      // register the values
      trans.EIRP = EIRP_dB;
      trans.G_T = G_T_rx_dB;
      trans.CN0 = CN0;
      trans.CN0_linear = pow(10, CN0 / 10.0);
    } else {
      // skip the link budget computation (for fixed noise case)
      trans.EIRP = 0.0;
      trans.G_T = 0.0;
      trans.CN0 = 0.0;
      trans.CN0_linear = 0.0;
    }

    trans.vis_occult = vis_occult;
    trans.vis_all = vis_all;

    return trans;
  }

  Real SpaceChannel::SolveLightTimeDelayRx(Transmitter *tx, Receiver *rx, Real t_rx) {
    Real tau = 0.0;
    auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx);
    auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);

    Real tau_prev = 0.0;
    Real rho = 0.0;
    int max_iter = 100;

    for (int n_iter = 0; n_iter < max_iter; n_iter++) {
      rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);
      rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();
      tau = rho / C;
      if (abs(tau - tau_prev) < 1e-12)
        break;
      else
        tau_prev = tau;
    }
    return tau;
  }

  Real SpaceChannel::SolveLightTimeDelayTx(Transmitter *tx, Receiver *rx, Real t_tx) {
    Real tau = 0.0;
    auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
    auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx);

    Real tau_prev = 0.0;
    Real rho = 0.0;
    int max_iter = 100;

    for (int n_iter = 0; n_iter < max_iter; n_iter++) {
      rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
      rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm();
      tau = rho / C;
      if (abs(tau - tau_prev) < 1e-12)
        break;
      else
        tau_prev = tau;
    }
    return tau;
  }

  Real SpaceChannel::ComputeFreeSpaceLossdB(Real dist, Real lambda) {
    Real path_loss = pow(4 * PI * dist / lambda, 2);
    Real path_loss_dB = 10 * log10(path_loss);
    return path_loss_dB;
  }

}  // namespace lupnt
