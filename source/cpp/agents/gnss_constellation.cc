/**
 * @file gnss_constellation.cpp
 * @author Stanford NAV LAB
 * @brief Gnss Constellation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/agents/gnss_constellation.h"

#include <filesystem>
#include <string>

#include "lupnt/core/constants.h"
#include "lupnt/core/file.h"
#include "lupnt/physics/orbit_state/tle.h"

namespace lupnt {

  void GnssConstellation::InitializeWithTle(GnssType gnss_type, std::string tle_filename,
                                            Ptr<IDynamics> dynamics, Ptr<GnssChannel> channel,
                                            double epoch0_tai) {
    SetChannel(channel);
    SetDynamics(dynamics);
    SetEpoch(epoch0_tai);
    LoadTleFile(gnss_type, tle_filename);
  }

  void GnssConstellation::LoadTleFile(GnssType gnss_type, std::string filename) {
    // std::filesystem::path path = GetFilePath(filename);

    for (auto tle : TLE::FromFile(filename)) {
      // std::cout << "Loaded satellite: " << tle.name << " PRN: " << tle.prn << std::endl;

      double sat_epoch = tle.epoch_tai;
      double dt_epoch = sat_epoch - epoch_;

      if (is_first_sat) {
        epoch_ = sat_epoch;
        is_first_sat = false;

        if (abs(dt_epoch) > 7 * SECS_DAY) {
          // warning message
          std::cout << "Warning: Satellite " << tle.prn << " has epoch " << sat_epoch
                    << " which is more than 7 days from the constellation epoch " << epoch_
                    << std::endl;
        }
      }

      double T = SECS_DAY / tle.mean_motion;

      // Classical orbital elements
      Real a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      Real e = tle.eccentricity;
      Real i = tle.inclination * RAD;
      Real Omega = tle.raan * RAD;
      Real w = tle.arg_perigee * RAD;
      Real rad_per_sec = tle.mean_motion * 2 * PI / SECS_DAY;  // TLE mean motion is in revs/day
      Real M = Wrap2Pi(tle.mean_anomaly * RAD + dt_epoch * rad_per_sec);  // [rad]

      // Convert to Cartesian
      ClassicalOE coe({a, e, i, Omega, w, M});
      coe.SetCoordSystem(Frame::GCRF);
      CartesianOrbitState cart = Classical2Cart(coe, GM_EARTH);
      auto state = MakePtr<CartesianOrbitState>(cart);

      // Create the spacecraft
      auto sat = MakePtr<Spacecraft>();
      sat->SetEpoch(sat_epoch);
      sat->SetDynamics(dynamics_);
      sat->SetOrbitState(state);
      sat->SetEpoch(epoch_);
      sat->SetBodyId(NaifId::EARTH);

      if (channel_) {
        auto transmitter = MakePtr<GnssTransmitter>(gnss_type, tle.prn);
        sat->AddDevice(transmitter);
        transmitter->SetAgent(sat);
        channel_->AddTransmitter(transmitter);
        transmitter->SetChannel(std::static_pointer_cast<SpaceChannel>(channel_));
      }

      satellites_.push_back(sat);
      gnss_types_.push_back(gnss_type);
    }

    // Propagate all satellites to the epoch
    for (auto sat : satellites_) {
      sat->Propagate(epoch_);
    }
  }

}  // namespace lupnt
