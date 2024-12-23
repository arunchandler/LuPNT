/**
 * @file link_measurement.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2024-04-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "radio_measurement.h"
#include "transmission.h"

namespace lupnt {

  enum class LinkMeasurementType { Range, RangeRate };

  struct LinkParams {
    // Receiver Parameters
    double freq;                 // carrier frequency [Hz]
    Modulation modulation_type;  // carrier type
    double B_L_chip;             // tracking loop noise bandwidth
    double Tc;                   // chip duration
    double B_L_carrier;          // carrier loop noise bandwidth
    double sigma_y_1s;           // one-sigma range noise [m]
    double m_R;                  // modulation index
    double T_I_doppler;          // Doppler integration time
    double T_I_range;            // range integration time (for open loop)
    double turnaround_ratio;     // Transponder turnaround ratio
    double elevation;            // Elevation angle [rad]

    // Agent Parameters
    Ptr<Agent> tx_agent;  // Transmitter (target) agent
    Ptr<Agent> rx_agent;  // Receiver agent

    double CN0_linear;  // Carrier-to-noise density [dB-Hz]
  };

  class LinkMeasurement {
  private:
    // state size
    int state_size_ow_ = 16;  // One way link state size (target rv + clock,  receiver rv + clock)
    int state_size_tw_ = 16;  // Two way link state size (target rv + receiver rv)

  protected:
    // True epochs
    Real epoch_tx_true_ = 0.0;  // True transmitted epoch
    Real epoch_rx_true_ = 0.0;  // True received epoch

    // Recorded epochs
    Real epoch_tx_recorded_ = 0.0;    // Recorded transmitted epoch
    Real epoch_rx_recorded_ = 0.0;    // Recorded received epoch
    Real epoch_rx_recorded_u_ = 0.0;  // Recorded received epoch (uplink, used for time sync)

    // Transmission
    ITransmission trans_ow_;                // One-way transmission
    std::vector<ITransmission> trans_tw_;   // Two-way transmission
    std::vector<ITransmission> trans_dow_;  // Dual one-way transmission

    // Parameters
    LinkParams linkparams_;

    // Flag if the link is generated
    bool one_way_generated_ = false;
    bool two_way_generated_ = false;
    bool dual_one_way_generated_ = false;

    // Random seed
    int seed_ = 0;

    // Occultation bodies
    std::vector<NaifId> occult_bodies_;
    VecXd occult_alt_;
    bool use_elev_mask_tx_ = false;
    bool use_elev_mask_rx_ = false;
    VecXd elev_masks_;

    // Use Fixed error for the measurements
    bool use_fixed_error_ = false;
    double range_sigma_fixed_ = 0.0;
    double range_rate_sigma_fixed_ = 0.0;

    // Visibility
    bool vis_ow_ = true;
    bool vis_tw_ = true;

    // bool
    bool use_open_loop_ = false;

    // Hardware delay (same for transmitter and receiver)
    Real hardware_delay_ = 1e-9;  // hardware delay [s]

  public:
    /**
     * @brief Construct a new Isl Measurement object using transmitters and
     * receivers
     *
     * @param occult_bodies  occulting bodies
     * @param occult_alt  occultation altitude
     * @param hardware_delay  hardware delay
     * @param link_type  link type  (e.g. "one-way", "two-way", "dual-one-way")
     */
    LinkMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt, VecXd elev_masks,
                    bool use_elev_mask_tx, bool use_elev_mask_rx, Real hardware_delay);

    /********************** Utils  *********************************/
    void SetLinkParams();
    inline void SetSeed(int seed) { seed_ = seed; }
    inline void SetFixedRangeError(double range_sigma) { range_sigma_fixed_ = range_sigma; }
    inline void SetFixedRangeRateError(double range_rate_sigma) {
      range_rate_sigma_fixed_ = range_rate_sigma;
    }
    inline void UseFixedError() { use_fixed_error_ = true; }
    inline void DisableFixedError() { use_fixed_error_ = false; }

    inline int GetSeed() const { return seed_; }
    inline LinkParams GetLinkParams() const { return linkparams_; }
    inline Real GetTxEpoch() const { return epoch_tx_true_; }
    inline Real GetRxEpoch() const { return epoch_rx_true_; }
    inline bool IsOneWayVisible() const { return vis_ow_; }
    inline bool IsTwoWayVisible() const { return vis_tw_; }
    inline Real GetRecordedEpochTx() const { return epoch_tx_recorded_; }
    inline Real GetRecordedEpochRx() const { return epoch_rx_recorded_; }
    inline Real GetRecordedEpochRxU() const { return epoch_rx_recorded_u_; }

    void Reset() {
      one_way_generated_ = false;
      two_way_generated_ = false;
      dual_one_way_generated_ = false;
      // reset transmission (one-way)
      trans_ow_ = ITransmission();
      trans_tw_.clear();
      trans_dow_.clear();

      // clear linkparams
      linkparams_ = LinkParams();

      // Reset recorded epochs
      epoch_tx_recorded_ = 0.0;
      epoch_rx_recorded_ = 0.0;
      epoch_rx_recorded_u_ = 0.0;
    }

    /********************** One way Link ***************************/

    /**
     * @brief Get the One Way Link object
     *
     * @param epoch_local     (local) epoch of the receiver or transmitter when the signal is
     * transmitted/received
     * @param tx              transmitter
     * @param rx              receiver
     * @param fixed_txrx      fixed transmitter or receiver time (tx or rx)
     * @return ITransmission
     */
    void GenerateOneWayLink(Real epoch_local, std::shared_ptr<Transmitter> &tx,
                            std::shared_ptr<Receiver> &rx, std::string fixed_txrx);

    void GenerateOneWayLinkAtRxEpoch(Real epoch, std::shared_ptr<Transmitter> &tx,
                                     std::shared_ptr<Receiver> &rx) {
      GenerateOneWayLink(epoch, tx, rx, "rx");
    };

    void GenerateOneWayLinkAtTxEpoch(Real epoch, std::shared_ptr<Transmitter> &tx,
                                     std::shared_ptr<Receiver> &rx) {
      GenerateOneWayLink(epoch, tx, rx, "tx");
    };

    VecX GetTrueOneWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types);

    /**
     * @brief Get the One Way Link Measurement
     *
     * @param epoch_rx_recorded  recorded epoch of the receiver (TAI)
     * @param epoch_ref          reference epoch
     * @param rv_tx              Transmitter position and velocity at the ref epoch
     * @param rv_rx              Receiver position and velocity at the ref epoch
     * @param clk_tx             Transmitter clock state at the ref epoch
     * @param clk_rx             Receiver clock state at the ref epoch
     * @param H_ow_rx            Jacobian matrix (n_meas x 16)
     * @param additional_delay   Additional delay for the rx time  (recorded = true + clk_bias +
     * additional_delay)
     * @param meas_types         measurement types
     * @param with_noise         (bool) add noise to the measurement
     * @param with_jacobian      (bool) compute jacobian
     * @return VecX
     */
    VecX GetOneWayLinkMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx, Vec6 rv_rx,
                                  Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx, Real additional_delay,
                                  std::vector<LinkMeasurementType> meas_types, bool with_noise,
                                  bool with_jacobian);

    /**
     * @brief Get the One Way Range Measurement
     *
     * @param epoch_rx_recorded  recorded epoch of the receiver (TAI)
     * @param epoch_ref          reference epoch
     * @param rv_tx              Transmitter position and velocity at the ref epoch
     * @param rv_rx              Receiver position and velocity at the ref epoch
     * @param clk_tx             Transmitter clock state at the ref epoch
     * @param clk_rx             Receiver clock state at the ref epoch
     * @param H_ow_rx            Jacobian matrix (1 x 16)
     * @param additional_delay   Additional delay for the rx time  (recorded = true + clk_bias +
     * additional_delay)
     * @param with_noise         (bool) add noise to the measurement
     * @param with_jacobian      (bool) compute jacobian
     * @return Real
     */
    Real GetOneWayRangeMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx, Vec6 rv_rx,
                                   Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx, Real additional_delay,
                                   bool with_noise, bool with_jacobian);

    /**
     * @brief Get the One Way Range Rate Measurement object
     *
     * @param epoch_rx_recorded  recorded epoch of the receiver (TAI)
     * @param epoch_ref          reference epoch
     * @param rv_tx              Transmitter position and velocity at the ref epoch
     * @param rv_rx              Receiver position and velocity at the ref epoch
     * @param clk_tx             Transmitter clock state at the ref epoch
     * @param clk_rx             Receiver clock state at the ref epoch
     * @param H_ow_rx            Jacobian matrix (1 x 16)
     * @param additional_delay   Additional delay for the rx time  (recorded = true + clk_bias +
     * additional_delay)
     * @param with_noise         (bool) add noise to the measurement
     * @param with_jacobian      (bool) compute jacobian
     * @return Real
     */
    Real GetOneWayRangeRateMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx,
                                       Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx,
                                       Real additional_delay, bool with_noise, bool with_jacobian);

    VecXd GetOneWayLinkNoise(std::vector<LinkMeasurementType> meas_types);
    double GetOneWayRangeNoise();
    double GetOneWayRangeRateNoise();

    /********************** Two way Link ***************************/

    /**
     * @brief Generate Two way link (receiver->target->receiver)
     *
     * @param epoch       epoch of the receiver or transmitter
     * @param tr_receiver    Transponder (receiver)
     * @param tr_target      Transponder (target)
     * @param txrx_fixed     fixed time of transmitter or receiver (tx or rx)
     * @return std::vector<ITransmission>
     */
    void GenerateTwoWayLink(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                            std::shared_ptr<Transponder> &tr_target, std::string txrx_fixed);

    void GenerateTwoWayLinkAtRxEpoch(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                                     std::shared_ptr<Transponder> &tr_target) {
      GenerateTwoWayLink(epoch, tr_receiver, tr_target, "rx");
    };

    void GenerateTwoWayLinkAtTxEpoch(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                                     std::shared_ptr<Transponder> &tr_target) {
      GenerateTwoWayLink(epoch, tr_receiver, tr_target, "tx");
    };

    /**
     * @brief Get the True Two Way Link Measurement object
     *
     * @param meas_types  vector of measurement types
     * @return VecX  measurement vector
     */
    VecX GetTrueTwoWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types);

    /**
     * @brief Get the Two Way Link Measurement (receiver->target->receiver)
     *
     * @param epoch_rx        Reception epoch (TAI) t_R
     * @param epoch_ref       Reference epoch
     * @param rv_receiver     Receiver position at ref epoch (w.r.t to
     * @param rv_target       Transmitter position at ref epoch (w.r.t to
     * @param clk_receiver    Receiver clock state at the ref epoch
     * @param clk_target      Transmitter clock state at the ref epoch
     * @param H_tw_rx         Jacobian matrix (n_meas x 16)
     * @param hardware_delay  Hardware delay at the target relay
     * @param meas_types      vector of measurement types
     * @param with_noise      (bool) add noise to the measurement
     * @param with_jacobian   (bool) compute jacobian
     * @return VecX
     */
    VecX GetTwoWayLinkMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver, Vec6 rv_target,
                                  Vec2 clk_receiver, Vec2 clk_target, MatXd &H_tw_rx,
                                  Real hardware_delay, std::vector<LinkMeasurementType> meas_types,
                                  bool with_noise, bool with_jacobian);

    /**
     * @brief Get the Two Way Range Measurement object
     *
     * @param epoch_rx        Reception epoch (TAI) t_R
     * @param epoch_ref       Reference epoch
     * @param rv_receiver     Receiver position at ref epoch (w.r.t to
     * @param rv_target       Transmitter position at ref epoch (w.r.t to
     * @param clk_receiver    Receiver clock state at the ref epoch
     * @param clk_target      Transmitter clock state at the ref epoch
     * @param H_tw_range      Jacobian matrix (1 x 16)
     * @param hardware_delay  Hardware delay at the target relay
     * @param with_noise      (bool) add noise to the measurement
     * @param with_jacobian   (bool) compute jacobian
     * @return Real
     */
    Real GetTwoWayRangeMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver, Vec6 rv_target,
                                   Vec2 clk_receiver, Vec2 clk_target, MatXd &H_tw_range,
                                   Real hardware_delay, bool with_noise, bool with_jacobian);

    /**
     * @brief Get the Two Way Range Rate Measurement object
     *
     * @param epoch_rx        Reception epoch (TAI) t_R
     * @param epoch_ref       Reference epoch
     * @param rv_receiver     Receiver position at ref epoch
     * @param rv_target       Transmitter position at ref epoch
     * @param clk_receiver    Receiver clock state at the ref epoch
     * @param clk_target      Transmitter clock state at the ref epoch
     * @param H_tw_rr         Jacobian matrix (1 x 16)
     * @param hardware_delay  Hardware delay at the target relay
     * @param with_noise      (bool) add noise to the measurement
     * @param with_jacobian   (bool) compute jacobian
     * @return Real
     */
    Real GetTwoWayRangeRateMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver,
                                       Vec6 rv_target, Vec2 clk_receiver, Vec2 clk_target,
                                       MatXd &H_tw_rr, Real hardware_delay, bool with_noise,
                                       bool with_jacobian);

    VecXd GetTwoWayLinkNoise(std::vector<LinkMeasurementType> meas_types);
    double GetTwoWayRangeNoise();
    double GetTwoWayRangeRateNoise();
  };

}  // namespace lupnt
