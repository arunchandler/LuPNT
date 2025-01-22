/**
 * @file groundstation_measurement.h
 * @author Stanford NAV Lab
 * @brief A set of measurement function for ground station observations
 * @version 0.1
 * @date 2024-12-03
 *
 * @copyright Copyright (c) 2024
 *
 */
#pragma once

#include "lupnt/agents/ground_station.h"
#include "lupnt/measurements/link_measurement.h"

namespace lupnt {

  struct DDOR_param {
    Real epoch_rx1_recorded = 0.0;
    Real epoch_rx2_recorded = 0.0;
    bool vis_ddor = false;
  };

  class GroundStationMeasurement : public LinkMeasurement {
  private:
    // Ground station ID
    int state_size_ow_gs_ = 8;  // One-way
    int state_size_tw_gs_ = 8;  // Two-way
    DDOR_param ddor_param_;

  public:
    GroundStationMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt, VecXd elev_masks,
                             Real hardware_delay)
        : LinkMeasurement(occult_bodies, occult_alt, elev_masks, true, true, hardware_delay) {
      use_open_loop_ = true;  // Open loop tracking for ground station
    }

    /********************** One way Down Link ***************************/
    Real GetOneWayRangeMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx, Vec6 rv_rx,
                                   Vec2 clk_tx, Vec2 clk_rx, Real additional_delay, bool with_noise,
                                   MatXd *H_ow_rx = nullptr);

    Real GetOneWayRangeRateMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx,
                                       Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx, Real additional_delay,
                                       bool with_noise, MatXd *H_ow_rx = nullptr);

    /********************** Two way Link ***************************/
    Real GetTwoWayRangeMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver, Vec6 rv_target,
                                   Vec2 clk_receiver, Vec2 clk_target, Real hardware_delay,
                                   bool with_noise, MatXd *H_tw_range = nullptr);

    Real GetTwoWayRangeRateMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver,
                                       Vec6 rv_target, Vec2 clk_receiver, Vec2 clk_target,
                                       Real hardware_delay, bool with_noise,
                                       MatXd *H_tw_rr = nullptr);

    /******** Delta DOR *****************/
    void GenerateDORLink(Real epoch_tx, Ptr<Transmitter> &tx, Ptr<Receiver> &rx1,
                         Ptr<Receiver> &rx2);

    Real GetDeltaDORMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_rx1, Vec6 rv_rx2, Vec6 rv_tx,
                                Vec2 clk_rx1, Vec2 clk_rx2, Vec2 clk_tx, Real hardware_delay,
                                bool with_noise, MatXd *H_dor = nullptr);
  };

}  // namespace lupnt
