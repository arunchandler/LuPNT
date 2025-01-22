/**
 * @file gnss_constellation.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/measurements/gnss_channel.h"
#include "lupnt/measurements/gnss_transmitter.h"
#include "lupnt/physics/orbit_state.h"
#include "spacecraft.h"

namespace lupnt {

  class GnssConstellation {
  private:
    std::vector<Ptr<Spacecraft>> satellites_;
    Ptr<IDynamics> dynamics_;
    Ptr<GnssChannel> channel_;
    Real epoch_;
    bool is_first_sat = true;
    std::vector<GnssType> gnss_types_;

  public:
    GnssConstellation() { is_first_sat = true; }

    void InitializeWithTle(GnssType gnss_type, std::string tle_filename, Ptr<IDynamics> dynamics,
                           Ptr<GnssChannel> channel, double epoch0_tai);

    void AddSatellitesWithTle(GnssType gnss_type, std::string tle_filename) {
      LoadTleFile(gnss_type, tle_filename);
    }

    void SetChannel(Ptr<GnssChannel> ch) { channel_ = ch; }

    void LoadTleFile(GnssType gnss_type, std::string filename);

    Real GetEpoch() { return epoch_; }

    void SetEpoch(double ep) {
      epoch_ = ep;
      for (auto sat : satellites_) sat->SetEpoch(ep);
    }

    void SetDynamics(Ptr<IDynamics> dyn) {
      dynamics_ = dyn;
      // set dynamics to all satellites
      for (auto sat : satellites_) sat->SetDynamics(dyn);
    }

    // Getters
    int GetNumSatellites() { return satellites_.size(); }
    int GetNumSatellites(GnssType gnss_type) {
      int count = 0;
      for (auto type : gnss_types_)
        if (type == gnss_type) count++;
      return count;
    }
    Ptr<Spacecraft> GetSatellite(int i) { return satellites_[i]; }
    Ptr<GnssChannel> GetChannel() { return channel_; }
    Ptr<IDynamics> GetDynamics() { return dynamics_; }

    // Methods
    void Propagate(Real epoch) {
      if (epoch == epoch_) return;
#pragma omp parallel for
      for (auto sat : satellites_) sat->Propagate(epoch);
      epoch_ = epoch;
    }
  };

}  // namespace lupnt
