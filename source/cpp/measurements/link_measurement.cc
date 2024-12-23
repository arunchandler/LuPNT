/**
 * @file link_measurement.cc
 * @author Stanford Nav Lab
 * @brief
 * @version 0.1
 * @date 2024-08-20
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/measurements/link_measurement.h"

#include "lupnt/core/constants.h"
#include "lupnt/measurements/radio_measurement.h"
#include "lupnt/measurements/space_channel.h"
#include "lupnt/physics/body.h"

namespace lupnt {

  LinkMeasurement::LinkMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt,
                                   VecXd elev_masks, bool use_elev_mask_tx, bool use_elev_mask_rx,
                                   Real hardware_delay)
      : occult_bodies_(occult_bodies),
        occult_alt_(occult_alt),
        elev_masks_(elev_masks),
        use_elev_mask_tx_(use_elev_mask_tx),
        use_elev_mask_rx_(use_elev_mask_rx),
        hardware_delay_(hardware_delay) {}

  void LinkMeasurement::SetLinkParams() {
    ITransmission trans, trans_u;
    // double G = 1.0;

    if (one_way_generated_) {
      trans = trans_ow_;

    } else if (two_way_generated_) {
      trans = trans_tw_[1];    // Downlink
      trans_u = trans_tw_[0];  // Uplink
    }

    // Link Parameters
    Ptr<Agent> tx_agent = trans.tx->GetAgent();
    Ptr<Agent> rx_agent = trans.rx->GetAgent();
    linkparams_.tx_agent = tx_agent;
    linkparams_.rx_agent = rx_agent;

    // Signal parameters
    ReceiverParam rx_param = trans.rx->rx_param_;
    linkparams_.freq = trans.tx->freq_tx;
    linkparams_.B_L_chip = rx_param.B_L_chip;
    linkparams_.B_L_carrier = rx_param.B_L_carrier;
    linkparams_.Tc = rx_param.modulation_type;
    linkparams_.T_I_doppler = rx_param.T_I_doppler;
    linkparams_.T_I_range = rx_param.T_I_range;
    linkparams_.elevation =

        // Singals
        linkparams_.CN0_linear = trans.CN0_linear;
  }

  /********************** One way Link ***************************/

  void LinkMeasurement::GenerateOneWayLink(Real epoch_local, std::shared_ptr<Transmitter> &tx,
                                           std::shared_ptr<Receiver> &rx, std::string txrx) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies_, occult_alt_);  // set occultation bodies
    sc.SetElevationMask(use_elev_mask_tx_, use_elev_mask_rx_, elev_masks_);  // set elevation mask
    Real t_tx_d, t_rx_d;

    bool compute_cn0 = true;
    if (use_fixed_error_) {
      compute_cn0 = false;
    }

    if (txrx == "rx") {
      t_rx_d = epoch_local
               - rx->GetAgent()->GetClockStateVecAtEpoch(epoch_local)(0);  // true received time
      trans_ow_ = sc.ComputeLinkBudget(tx, rx, t_rx_d, "rx", compute_cn0);
      t_tx_d = trans_ow_.t_tx;  // true transmission time
      epoch_tx_recorded_ = t_tx_d + tx->GetAgent()->GetClockStateVecAtEpoch(t_tx_d)(0);
    } else if (txrx == "tx") {
      t_tx_d = epoch_local
               - tx->GetAgent()->GetClockStateVecAtEpoch(epoch_local)(0);  // true transmission time
      trans_ow_ = sc.ComputeLinkBudget(tx, rx, t_tx_d, "tx", compute_cn0);
      t_rx_d = trans_ow_.t_rx;  // true received time
      epoch_rx_recorded_
          = rx->GetAgent()->GetClockStateVecAtEpoch(t_rx_d)(0);  // recorded received time
    }
    one_way_generated_ = true;

    vis_ow_ = trans_ow_.vis_all;  // set visibility

    // Set time
    epoch_tx_true_ = t_tx_d;
    epoch_rx_true_ = t_rx_d;

    // Link Parameters
    SetLinkParams();
  }

  VecX LinkMeasurement::GetTrueOneWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types) {
    if (!one_way_generated_) {
      std::cerr << "Error: One way link not generated" << std::endl;
      return VecX::Zero(0);
    }
    if (meas_types.size() == 0) {
      std::cerr << "Error: No measurement types provided" << std::endl;
      return VecX::Zero(0);
    }
    if (vis_ow_ == false) {  // No visibility
      return VecX::Zero(0);
    }

    VecX rv_tx = linkparams_.tx_agent->GetRvStateAtEpoch(epoch_rx_true_);
    VecX rv_rx = linkparams_.rx_agent->GetRvStateAtEpoch(epoch_rx_true_);
    VecX clk_tx = linkparams_.tx_agent->GetClockStateVecAtEpoch(epoch_rx_true_);
    VecX clk_rx = linkparams_.rx_agent->GetClockStateVecAtEpoch(epoch_rx_true_);

    MatXd H_ow_rx(2, 8);  // temporary, won't be used

    VecX y = GetOneWayLinkMeasurement(epoch_rx_recorded_, epoch_rx_true_, rv_tx, rv_rx, clk_tx,
                                      clk_rx, H_ow_rx, 0.0, meas_types, true, false);

    return y;
  }

  VecX LinkMeasurement::GetOneWayLinkMeasurement(Real epoch_rx_recorded, Real epoch_ref, Vec6 rv_tx,
                                                 Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                                 MatXd &H_ow_rx, Real additional_delay,
                                                 std::vector<LinkMeasurementType> meas_types,
                                                 bool with_noise, bool with_jacobian) {
    Real rho_ow, rho_ow_rate;
    // double sigma_ow = 0.0;
    // double sigma_ow_rate = 0.0;

    if (vis_ow_ == false) {  // No visibility
      return VecX::Zero(0);
    }

    int meas_size = meas_types.size();
    int state_size = state_size_ow_;

    VecX z(meas_size);
    H_ow_rx.resize(meas_size, state_size);
    MatXd H_ow_range = MatXd::Zero(1, state_size);
    MatXd H_ow_rangerate = MatXd::Zero(1, state_size);

    int idx = 0;

    if (meas_size == 0) {
      return VecX::Zero(0);
    }

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case LinkMeasurementType::Range: {
          rho_ow = GetOneWayRangeMeasurement(epoch_rx_recorded, epoch_ref, rv_tx, rv_rx, clk_tx,
                                             clk_rx, H_ow_range, additional_delay, with_noise,
                                             with_jacobian);
          z(idx) = rho_ow;
          H_ow_rx.row(idx) = H_ow_range;
          idx++;
        }
        case LinkMeasurementType::RangeRate: {
          rho_ow_rate = GetOneWayRangeRateMeasurement(epoch_rx_recorded, epoch_ref, rv_tx, rv_rx,
                                                      clk_tx, clk_rx, H_ow_rangerate,
                                                      additional_delay, with_noise, with_jacobian);
          z(idx) = rho_ow_rate;
          H_ow_rx.row(idx) = H_ow_rangerate;
          idx++;
        }
        default: break;
      }
    }

    return z;
  }

  Real LinkMeasurement::GetOneWayRangeMeasurement(Real epoch_rx_recorded, Real epoch_ref,
                                                  Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                                  MatXd &H_ow_rx, Real additional_delay,
                                                  bool with_noise, bool with_jacobian) {
    // return

    auto func = [epoch_rx_recorded, epoch_ref, additional_delay, this](
                    const Vec6 rv_tx_in, const Vec2 clk_tx_in, const Vec6 rv_rx_in,
                    const Vec2 clk_rx_in) {
      Real owr = ComputeOneWayRangeLTR(epoch_rx_recorded, epoch_ref, rv_tx_in, rv_rx_in, clk_tx_in,
                                       clk_rx_in, linkparams_.tx_agent, linkparams_.rx_agent,
                                       additional_delay, false);
      return owr;
    };

    // Without lighttime delay
    // auto func = [rv_tx, clk_tx](const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
    //   Real owr = ComputePseudorange(
    //       rv_tx.head(3), rv_rx_in.head(3), clk_tx(0), clk_rx_in(0), 0.0);
    //   return owr;
    // };

    if (vis_ow_ == false) {  // No visibility
      return 0.0;
    }

    // break the computational graph relations before taking the jacobian
    Vec6 rv_tx_tmp = rv_tx.cast<double>();
    Vec2 clk_tx_tmp = clk_tx.cast<double>();
    Vec6 rv_rx_tmp = rv_rx.cast<double>();
    Vec2 clk_rx_tmp = clk_rx.cast<double>();
    Real rho_ow = 0.0;

    if (with_jacobian) {
      VecXd H_ow_vec = gradient(func, wrt(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp),
                                at(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp), rho_ow);
      // Convert to (1, 8) matrix
      H_ow_rx.row(0) = H_ow_vec.transpose();
    } else {
      rho_ow = func(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp);
    }

    if (with_noise) {
      double sigma_ow = GetOneWayRangeNoise();
      rho_ow += SampleRandNormal(0.0, sigma_ow, seed_);
    }

    return rho_ow;
  }

  Real LinkMeasurement::GetOneWayRangeRateMeasurement(Real epoch_rx_recorded, Real epoch_ref,
                                                      Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                                      Vec2 clk_rx, MatXd &H_ow_rx,
                                                      Real hardware_delay, bool with_noise,
                                                      bool with_jacobian) {
    auto func = [epoch_rx_recorded, epoch_ref, rv_tx, clk_tx, hardware_delay, this](
                    const Vec6 rv_tx_in, const Vec2 clk_tx_in, const Vec6 rv_rx_in,
                    const Vec2 clk_rx_in) {
      Real owrr = ComputeOneWayRangeRateLTR(
          epoch_rx_recorded, epoch_ref, rv_tx_in, rv_rx_in, clk_tx_in, clk_rx_in,
          linkparams_.tx_agent, linkparams_.rx_agent, hardware_delay, linkparams_.T_I_doppler);
      return owrr;
    };

    if (vis_ow_ == false) {  // No visibility
      return 0.0;
    }

    // break the computational graph relations before taking the jacobian
    Vec6 rv_tx_tmp = rv_rx.cast<double>();
    Vec2 clk_tx_tmp = clk_rx.cast<double>();
    Vec6 rv_rx_tmp = rv_rx.cast<double>();
    Vec2 clk_rx_tmp = clk_rx.cast<double>();
    Real rho_ow_rate = 0.0;

    if (with_jacobian) {
      VecXd H_ow_rate_vec = gradient(func, wrt(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp),
                                     at(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp), rho_ow_rate);

      // Convert to (1, nstate) matrix
      H_ow_rx.row(0) = H_ow_rate_vec.transpose();

    } else {
      rho_ow_rate = func(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp);
    }

    if (with_noise) {
      double sigma_ow_rate = GetOneWayRangeRateNoise();
      rho_ow_rate += SampleRandNormal(0.0, sigma_ow_rate, seed_);
    }

    return rho_ow_rate;
  }

  /* Noise models for One-way */
  VecXd LinkMeasurement::GetOneWayLinkNoise(std::vector<LinkMeasurementType> meas_types) {
    VecXd noise_std_vec(meas_types.size());
    int idx = 0;

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case LinkMeasurementType::Range: {
          noise_std_vec(idx) = GetOneWayRangeNoise();
          idx += 1;
          break;
        }
        case LinkMeasurementType::RangeRate: {
          noise_std_vec(idx) = GetOneWayRangeRateNoise();
          idx += 1;
          break;
        }
        default: break;
      }
    }

    return noise_std_vec;
  }

  double LinkMeasurement::GetOneWayRangeNoise() {
    double sigma_ow = 0.0;

    if (use_fixed_error_) {
      sigma_ow = range_sigma_fixed_;
    } else {
      // For one-way link, double the range error of the two way link
      if (!use_open_loop_) {
        sigma_ow = 2
                   * ComputePnRangeErrorCTL(linkparams_.CN0_linear, linkparams_.B_L_chip,
                                            linkparams_.Tc, linkparams_.modulation_type);
      } else {
        sigma_ow = 2
                   * ComputePnRangeErrorOL(linkparams_.CN0_linear, linkparams_.T_I_range,
                                           linkparams_.Tc, linkparams_.modulation_type);
      }
    }
    return sigma_ow;
  }

  double LinkMeasurement::GetOneWayRangeRateNoise() {
    double sigma_ow_rate = 0.0;

    if (use_fixed_error_) {
      sigma_ow_rate = range_rate_sigma_fixed_;
    } else {
      sigma_ow_rate = ComputeRangeRateErrorOneWay(linkparams_.B_L_carrier, linkparams_.freq,
                                                  linkparams_.Tc, linkparams_.T_I_doppler,
                                                  trans_ow_.CN0_linear, linkparams_.sigma_y_1s,
                                                  linkparams_.modulation_type, linkparams_.m_R);
    }
    return sigma_ow_rate;
  }

  /********************** Two way Link ***************************/

  void LinkMeasurement::GenerateTwoWayLink(Real epoch_local,
                                           std::shared_ptr<Transponder> &tr_receiver,
                                           std::shared_ptr<Transponder> &tr_target,
                                           std::string txrx) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies_, occult_alt_);  // set occultation bodies
    sc.SetElevationMask(use_elev_mask_tx_, use_elev_mask_rx_, elev_masks_);  // set elevation mask

    Real t_tx_u, t_rx_u, t_tx_d, t_rx_d;
    ITransmission trans_d, trans_u;

    // tr_target -> tr_receiver (downlink)
    std::shared_ptr<Transmitter> tx_d = tr_target->GetTransmitter();
    std::shared_ptr<Receiver> rx_d = tr_receiver->GetReceiver();

    // tr_receiver -> tr_target (uplink)
    std::shared_ptr<Transmitter> tx_u = tr_receiver->GetTransmitter();
    std::shared_ptr<Receiver> rx_u = tr_target->GetReceiver();

    // Get clock offsets

    bool compute_cn0 = true;
    if (use_fixed_error_) {
      compute_cn0 = false;
    }

    if (txrx == "rx") {
      Real delay_rx_d = tr_receiver->GetAgent()->GetClockStateVecAtEpoch(epoch_local)(0);
      t_rx_d = epoch_local - delay_rx_d;  // true received time
      trans_d = sc.ComputeLinkBudget(tx_d, rx_d, t_rx_d, "rx", compute_cn0);
      t_tx_d = trans_d.t_tx;
      t_rx_u = t_tx_d - hardware_delay_;  // true received time
      trans_u = sc.ComputeLinkBudget(tx_u, rx_u, t_rx_u, "rx", compute_cn0);
      t_tx_u = trans_u.t_tx;  // true transmission time

      // recorded time
      epoch_tx_recorded_ = t_tx_u + tr_receiver->GetAgent()->GetClockStateVecAtEpoch(t_tx_u)(0);
      epoch_rx_recorded_ = epoch_local;
      epoch_rx_recorded_u_ = t_rx_u + tr_target->GetAgent()->GetClockStateVecAtEpoch(t_rx_u)(0);

    } else if (txrx == "tx") {
      Real delay_tx_u = tr_receiver->GetAgent()->GetClockStateVecAtEpoch(epoch_local)(0);
      t_tx_u = epoch_local - delay_tx_u;  // true transmission
      trans_u = sc.ComputeLinkBudget(tx_u, rx_u, t_tx_u, "tx", compute_cn0);
      t_rx_u = trans_d.t_rx;
      t_tx_d = t_rx_u - hardware_delay_;
      trans_d = sc.ComputeLinkBudget(tx_d, rx_d, t_tx_d, "tx", compute_cn0);
      t_rx_d = trans_d.t_rx;

      // recoreded time
      epoch_tx_recorded_ = epoch_local;
      epoch_rx_recorded_ = t_rx_d + tr_receiver->GetAgent()->GetClockStateVecAtEpoch(t_rx_d)(0);
      epoch_rx_recorded_u_ = t_rx_u + tr_target->GetAgent()->GetClockStateVecAtEpoch(t_rx_u)(0);
    }

    two_way_generated_ = true;

    // Set visibility
    vis_tw_ = trans_u.vis_all & trans_d.vis_all;

    // Set time
    epoch_tx_true_ = t_tx_u;
    epoch_rx_true_ = t_rx_d;

    // Set link
    trans_tw_.push_back(trans_u);
    trans_tw_.push_back(trans_d);

    // Link Parameters -------------------------------------------------
    linkparams_.turnaround_ratio = tr_target->turnaround_ratio;

    SetLinkParams();
  }

  VecX LinkMeasurement::GetTrueTwoWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types) {
    if (!two_way_generated_) {
      std::cerr << "Error: Two way link not generated" << std::endl;
      return VecX::Zero(0);
    }

    if (meas_types.size() == 0) {
      std::cerr << "Error: No measurement types provided" << std::endl;
      return VecX::Zero(0);
    }

    if (vis_tw_ == false) {  // No visibility
      return VecX::Zero(0);
    }

    VecX rv_receiver = linkparams_.rx_agent->GetRvStateAtEpoch(epoch_rx_true_);
    VecX rv_target = linkparams_.tx_agent->GetRvStateAtEpoch(epoch_rx_true_);
    Vec2 clk_receiver = linkparams_.rx_agent->GetClockStateVecAtEpoch(epoch_rx_true_);
    Vec2 clk_target = linkparams_.tx_agent->GetClockStateVecAtEpoch(epoch_rx_true_);

    MatXd H_tw_rx(2, state_size_tw_);  // temporary, won't be used

    VecX y = GetTwoWayLinkMeasurement(epoch_rx_recorded_, epoch_rx_true_, rv_receiver, rv_target,
                                      clk_receiver, clk_target, H_tw_rx, hardware_delay_,
                                      meas_types, true, false);

    return y;
  }

  VecX LinkMeasurement::GetTwoWayLinkMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver,
                                                 Vec6 rv_target, Vec2 clk_receiver, Vec2 clk_target,
                                                 MatXd &H_tw_rx, Real hardware_delay,
                                                 std::vector<LinkMeasurementType> meas_types,
                                                 bool with_noise, bool with_jacobian) {
    // double sigma_tw = 0.0;
    // double sigma_tw_rate = 0.0;

    int meas_size = meas_types.size();
    int state_size = state_size_tw_;

    VecX z(meas_size);
    H_tw_rx.resize(meas_size, state_size);

    MatXd H_tw_range = MatXd::Zero(1, state_size);
    MatXd H_tw_rangerate = MatXd::Zero(1, state_size);

    int idx = 0;

    if (meas_size == 0) {
      return VecX::Zero(0);
    }
    if (vis_tw_ == false) {  // No visibility
      return VecX::Zero(0);
    }

    for (int i = 0; i < meas_size; i++) {
      LinkMeasurementType meas_type = meas_types[i];

      switch (meas_type) {
        case LinkMeasurementType::Range: {
          Real rho_tw = GetTwoWayRangeMeasurement(epoch_rx, epoch_ref, rv_receiver, rv_target,
                                                  clk_receiver, clk_target, H_tw_range,
                                                  hardware_delay, with_noise, with_jacobian);
          z(idx) = rho_tw;
          H_tw_rx.row(idx) = H_tw_range;
          idx++;
          break;
        }
        case LinkMeasurementType::RangeRate: {
          Real rho_tw_rate = GetTwoWayRangeRateMeasurement(
              epoch_rx, epoch_ref, rv_receiver, rv_target, clk_receiver, clk_target, H_tw_rangerate,
              hardware_delay, with_noise, with_jacobian);
          z(idx) = rho_tw_rate;
          H_tw_rx.row(idx) = H_tw_rangerate;
          idx++;
          break;
        }
        default: break;
      }
    }

    return z;
  }

  Real LinkMeasurement::GetTwoWayRangeMeasurement(Real epoch_rx, Real epoch_ref, Vec6 rv_receiver,
                                                  Vec6 rv_target, Vec2 clk_receiver,
                                                  Vec2 clk_target, MatXd &H_tw_range,
                                                  Real hardware_delay, bool with_noise,
                                                  bool with_jacobian) {
    // return

    auto func = [epoch_ref, epoch_rx, hardware_delay, this](
                    const Vec6 rv_target_in, const Vec2 clk_target_in, const Vec6 rv_receiver_in,
                    const Vec2 clk_receiver_in) {
      Real twr = ComputeTwoWayRangeLTR(epoch_rx, epoch_ref, rv_target_in, rv_receiver_in,
                                       clk_target_in, clk_receiver_in, linkparams_.tx_agent,
                                       linkparams_.rx_agent, hardware_delay, 0.0);
      return twr;
    };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_target_tmp = rv_target.cast<double>();
    Vec6 rv_receiver_tmp = rv_receiver.cast<double>();
    Vec2 clk_target_tmp = clk_target.cast<double>();
    Vec2 clk_receiver_tmp = clk_receiver.cast<double>();
    Real rho_tw = 0.0;

    if (with_jacobian) {
      VecXd H_tw_vec
          = gradient(func, wrt(rv_target_tmp, clk_target_tmp, rv_receiver_tmp, clk_receiver_tmp),
                     at(rv_target_tmp, clk_target_tmp, rv_receiver_tmp, clk_receiver_tmp), rho_tw);

      // Convert to (1, state) matrix
      H_tw_range.row(0) = H_tw_vec.transpose();
    } else {
      rho_tw = func(rv_target_tmp, clk_target_tmp, rv_receiver_tmp, clk_receiver_tmp);
    }

    if (with_noise) {
      double sigma_tw = GetTwoWayRangeNoise();
      double noise_tw = SampleRandNormal(0.0, sigma_tw, seed_);
      rho_tw += noise_tw;

      // // debug message
      // std::cout << " <Range Computation>" << std::endl;
      // std::cout << " Target Position    : " << rv_target_tmp.transpose() << std::endl;
      // std::cout << " Receiver Position  : " << rv_receiver_tmp.transpose() << std::endl;
      // std::cout << " Inst Distance      : " << (rv_target_tmp.head(3) -
      // rv_receiver_tmp.head(3)).norm() << std::endl; std::cout << " Light-Time Distance: " <<
      // rho_tw << std::endl; std::cout << " Noise              : " << noise_tw << std::endl;
      // std::cout << " " << std::endl;
    }
    return rho_tw;
  }

  Real LinkMeasurement::GetTwoWayRangeRateMeasurement(Real epoch_rx, Real epoch_ref,
                                                      Vec6 rv_receiver, Vec6 rv_target,
                                                      Vec2 clk_receiver, Vec2 clk_target,
                                                      MatXd &H_tw_rr, Real hardware_delay,
                                                      bool with_noise, bool with_jacobian) {
    // Function to compute the two way range rate
    auto func = [epoch_rx, epoch_ref, rv_target, hardware_delay, this](
                    const Vec6 rv_target_in, const Vec2 clk_target_in, const Vec6 rv_receiver_in,
                    const Vec2 clk_receiver_in) {
      Real owrr = ComputeTwoWayRangeRateLTR(
          epoch_rx, epoch_ref, rv_target_in, rv_receiver_in, clk_target_in, clk_receiver_in,
          linkparams_.tx_agent, linkparams_.rx_agent, hardware_delay, linkparams_.T_I_doppler);
      return owrr;
    };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_target_tmp = rv_target.cast<double>();
    Vec2 clk_target_tmp = clk_target.cast<double>();
    Vec6 rv_rx_tmp = rv_receiver.cast<double>();
    Vec2 clk_rx_tmp = clk_receiver.cast<double>();
    Real rho_tw_rate = 0.0;

    if (with_jacobian) {
      VecXd H_tw_rate_vec
          = gradient(func, wrt(rv_target_tmp, clk_target_tmp, rv_rx_tmp, clk_rx_tmp),
                     at(rv_target_tmp, clk_target_tmp, rv_rx_tmp, clk_rx_tmp), rho_tw_rate);

      // Convert to (1, 8) matrix
      H_tw_rr.row(0) = H_tw_rate_vec.transpose();

    } else {
      rho_tw_rate = func(rv_target_tmp, clk_target_tmp, rv_rx_tmp, clk_rx_tmp);
    }

    if (with_noise) {
      double sigma_tw_rate = GetTwoWayRangeRateNoise();
      double noise_tw_rate = SampleRandNormal(0.0, sigma_tw_rate, seed_);
      rho_tw_rate += noise_tw_rate;
    }

    return rho_tw_rate;
  }

  /* Noise models for Two-way */
  VecXd LinkMeasurement::GetTwoWayLinkNoise(std::vector<LinkMeasurementType> meas_types) {
    VecXd noise_std_vec(meas_types.size());
    int idx = 0;

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case LinkMeasurementType::Range: {
          noise_std_vec(idx) = GetTwoWayRangeNoise();
          idx += 1;
          break;
        }
        case LinkMeasurementType::RangeRate: {
          noise_std_vec(idx) = GetTwoWayRangeRateNoise();
          idx += 1;
          break;
        }
        default: break;
      }
    }

    return noise_std_vec;
  }

  double LinkMeasurement::GetTwoWayRangeNoise() {
    double sigma_tw = 0.0;

    if (use_fixed_error_) {
      sigma_tw = range_sigma_fixed_;
    } else {
      if (!use_open_loop_) {
        sigma_tw = ComputePnRangeErrorCTL(linkparams_.CN0_linear, linkparams_.B_L_chip,
                                          linkparams_.Tc, linkparams_.modulation_type);
      } else {
        sigma_tw = ComputePnRangeErrorOL(linkparams_.CN0_linear, linkparams_.T_I_range,
                                         linkparams_.Tc, linkparams_.modulation_type);
      }
    }

    return sigma_tw;
  }

  double LinkMeasurement::GetTwoWayRangeRateNoise() {
    double sigma_tw_rate = 0.0;

    if (use_fixed_error_) {
      sigma_tw_rate = range_rate_sigma_fixed_;
    } else {
      sigma_tw_rate = ComputeRangeRateErrorTwoWay(
          linkparams_.B_L_carrier, linkparams_.freq, linkparams_.Tc, linkparams_.T_I_doppler,
          linkparams_.CN0_linear, linkparams_.sigma_y_1s, linkparams_.turnaround_ratio,
          linkparams_.modulation_type, linkparams_.m_R);
    }
    return sigma_tw_rate;
  }

}  // namespace lupnt
