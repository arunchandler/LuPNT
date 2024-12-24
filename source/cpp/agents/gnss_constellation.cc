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
#include "lupnt/core/user_file_path.h"
#include "lupnt/physics/orbit_state/tle.h"

namespace lupnt {

  void GnssConstellation::InitializeWithTle(std::string gnss_type, std::string tle_filename, Ptr<IDynamics> dynamics, Ptr<GnssChannel> channel) {
    SetChannel(channel);
    SetDynamics(dynamics);
    LoadTleFile(gnss_type, tle_filename);
  }

  void GnssConstellation::LoadTleFile(std::string_view gnss_type, std::string filename) {
    // std::filesystem::path path = GetFilePath(filename);

    bool is_first = true;

    for (auto tle : TLE::FromFile(filename)) {
      // extract epoch
      if (is_first) {
        epoch_ = tle.epoch_tai;
      }
      // std::cout << "Loaded satellite: " << tle.name << " PRN: " << tle.prn << std::endl;

      double dt_epoch = epoch_ - tle.epoch_tai;

      double T = SECS_DAY / tle.mean_motion;

      // Classical orbital elements
      Real a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      Real e = tle.eccentricity;
      Real i = tle.inclination * RAD;
      Real Omega = tle.raan * RAD;
      Real w = tle.arg_perigee * RAD;
      Real rad_per_sec = tle.mean_motion * 2 * PI / SECS_DAY;  // TLE mean motion is in revs/day
      Real M = Wrap2Pi(tle.mean_anomaly * RAD + dt_epoch * rad_per_sec);  // [rad]

      ClassicalOE coe({a, e, i, Omega, w, M});
      coe.SetCoordSystem(Frame::GCRF);
      CartesianOrbitState cart = Classical2Cart(coe, GM_EARTH);
      auto state = MakePtr<CartesianOrbitState>(cart);

      // Create the spacecraft
      auto sat = MakePtr<Spacecraft>();
      sat->SetDynamics(dynamics_);
      sat->SetOrbitState(state);
      sat->SetEpoch(epoch_);
      sat->SetBodyId(NaifId::EARTH);

      if (channel_) {
        auto transmitter = MakePtr<GnssTransmitter>(std::string(gnss_type), tle.prn);
        sat->AddDevice(transmitter);
        transmitter->SetAgent(sat);
        channel_->AddTransmitter(transmitter);
        transmitter->SetChannel(std::static_pointer_cast<SpaceChannel>(channel_));
      }

      satellites_.push_back(sat);

      is_first = false;
    }
  }

}  // namespace lupnt
