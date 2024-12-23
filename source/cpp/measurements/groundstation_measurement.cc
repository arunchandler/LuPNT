/**
 * @file groundstation_measurement.cc
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2024-12-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/measurements/groundstation_measurement.h"

namespace lupnt {

  Real GroundStationMeasurement::GetOneWayRangeMeasurement(Real epoch_rx_recorded, Real epoch_ref,
                                                           Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                                           Vec2 clk_rx, MatXd &H_ow_rx,
                                                           Real additional_delay, bool with_noise,
                                                           bool with_jacobian) {
    MatXd H_ow_rx_tmp = MatXd::Zero(1, 16);
    LinkMeasurement::GetOneWayRangeMeasurement(epoch_rx_recorded, epoch_ref, rv_tx, rv_rx, clk_tx,
                                               clk_rx, H_ow_rx_tmp, additional_delay, with_noise,
                                               with_jacobian);

    // extract the first 8 states (tx rv + clk)
    H_ow_rx.resize(1, 8);
    H_ow_rx = H_ow_rx_tmp.block(0, 0, 1, 8);

    return H_ow_rx_tmp(0, 8);
  }

  Real GroundStationMeasurement::GetOneWayRangeRateMeasurement(
      Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
      MatXd &H_ow_rx, Real additional_delay, bool with_noise, bool with_jacobian) {
    MatXd H_ow_rx_tmp = MatXd::Zero(1, 16);
    LinkMeasurement::GetOneWayRangeRateMeasurement(epoch_rx_recorded, epoch_ref, rv_tx, rv_rx,
                                                   clk_tx, clk_rx, H_ow_rx_tmp, additional_delay,
                                                   with_noise, with_jacobian);

    // extract the first 8 states (tx rv + clk)
    H_ow_rx.resize(1, 8);
    H_ow_rx = H_ow_rx_tmp.block(0, 0, 1, 8);

    return H_ow_rx_tmp(0, 8);
  }

  Real GroundStationMeasurement::GetTwoWayRangeMeasurement(Real epoch_rx, Real epoch_ref,
                                                           Vec6 rv_receiver, Vec6 rv_target,
                                                           Vec2 clk_receiver, Vec2 clk_target,
                                                           MatXd &H_tw_range, Real hardware_delay,
                                                           bool with_noise, bool with_jacobian) {
    MatXd H_tw_rx_tmp = MatXd::Zero(1, 16);
    LinkMeasurement::GetTwoWayRangeMeasurement(epoch_rx, epoch_ref, rv_receiver, rv_target,
                                               clk_receiver, clk_target, H_tw_rx_tmp,
                                               hardware_delay, with_noise, with_jacobian);

    // extract the first 8 states (tx rv + clk)
    H_tw_range.resize(1, 8);
    H_tw_range = H_tw_rx_tmp.block(0, 0, 1, 8);

    return H_tw_rx_tmp(0, 8);
  }

  Real GroundStationMeasurement::GetTwoWayRangeRateMeasurement(
      Real epoch_rx, Real epoch_ref, Vec6 rv_receiver, Vec6 rv_target, Vec2 clk_receiver,
      Vec2 clk_target, MatXd &H_tw_rr, Real hardware_delay, bool with_noise, bool with_jacobian) {
    MatXd H_tw_rx_tmp = MatXd::Zero(1, 16);
    LinkMeasurement::GetTwoWayRangeRateMeasurement(epoch_rx, epoch_ref, rv_receiver, rv_target,
                                                   clk_receiver, clk_target, H_tw_rx_tmp,
                                                   hardware_delay, with_noise, with_jacobian);

    // extract the first 8 states (tx rv + clk)
    H_tw_rr.resize(1, 8);
    H_tw_rr = H_tw_rx_tmp.block(0, 0, 1, 8);

    return H_tw_rx_tmp(0, 8);
  }

  void GroundStationMeasurement::GenerateDORLink(Real epoch_tx, std::shared_ptr<Transmitter> &tx,
                                                 std::shared_ptr<Receiver> &rx1,
                                                 std::shared_ptr<Receiver> &rx2) {
    // Generate the one way link
    Real tx_clk_bias = tx->GetAgent()->GetClockStateVecAtEpoch(epoch_tx)(0);
    Real epoch_tx_local = epoch_tx + tx_clk_bias;
    LinkMeasurement::GenerateOneWayLink(epoch_tx_local, tx, rx1, "tx");
    bool vis_ow1 = vis_ow_;
    Real epoch_rx1_rec = epoch_rx_recorded_;

    // Generate the second one way link
    LinkMeasurement::GenerateOneWayLink(epoch_tx_local, tx, rx2, "tx");
    bool vis_ow2 = vis_ow_;
    Real epoch_rx2_rec = epoch_rx_recorded_;

    // store data
    ddor_param_.epoch_rx1_recorded = epoch_rx1_rec;
    ddor_param_.epoch_rx2_recorded = epoch_rx2_rec;
    ddor_param_.vis_ddor = vis_ow1 && vis_ow2;
  }

  Real GroundStationMeasurement::GetDeltaDORMeasurement(Real epoch_rx1, Real epoch_ref, Vec6 rv_rx1,
                                                        Vec6 rv_rx2, Vec6 rv_tx, Vec2 clk_rx1,
                                                        Vec2 clk_rx2, Vec2 clk_tx, MatXd &H_dor,
                                                        Real hardware_delay, bool with_noise,
                                                        bool with_jacobian) {}

}  // namespace lupnt
