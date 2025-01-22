/**
 * @file ground_stations.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/coordinates.h"

namespace lupnt {

  struct GroundStationData {
    int gs_id;
    std::string name;
    std::string location;
    double dish_size;
    double latitude;   // deg
    double longitude;  // deg
    double altitude_m;
    double eirp_min;  // dBW
    double eirp_max;
    double uplink_freq_min;  // MHz
    double uplink_freq_max;
    double downlink_freq_min;
    double downlink_freq_max;
    double gt_dBk;  // dB/K
    bool angle_meas;
    bool doppler_meas;
    bool range_meas;
    std::string freq_band;  // S, X, K, Ka
  };

  typedef std::map<int, GroundStationData> GroundStationDataMap;

  GroundStationDataMap LoadGroundStationData();

  class GroundStation : public Agent {
  protected:
    NaifId body_id_ = NaifId::EARTH;
    Real latitude_ = 0.0;   // [rad]
    Real longitude_ = 0.0;  // [rad]
    Real altitude_ = 0.0;   // [km]
    Frame frame_ = Frame::NONE;

  public:
    GroundStation(GroundStationData gs_data) : Agent() {
      name_ = gs_data.name;
      latitude_ = gs_data.latitude * RAD;
      longitude_ = gs_data.longitude * RAD;
      altitude_ = gs_data.altitude_m * 1e-3;

      SetIsBodyFixed(true);
      SetBodyId(body_id_);
      // Create body for corresponding body_id
      BodyData body_data = GetBodyData(body_id_);

      SetDynamics(MakePtr<SurfaceStaticDynamics>(body_id_, body_data.fixed_frame));

      // convert (lat, lon) to ECEF
      Vec3 lla = {latitude_, longitude_, altitude_};
      Vec3 pos_f = LatLonAlt2Cart(lla, body_data.R, body_data.flattening);

      // Set Cartesian State
      Frame frame = body_data.fixed_frame;
      frame_ = frame;

      Vec6 rv = {pos_f(0), pos_f(1), pos_f(2), 0.0, 0.0, 0.0};
      SetRvState(MakePtr<CartesianOrbitState>(rv, frame));
      SetDynamicsFrame(frame);
    };

    double GetLatitudeDouble() { return latitude_.val(); }
    double GetLongitudeDouble() { return longitude_.val(); }
    double GetAltitudeDouble() { return altitude_.val(); }
    VecXd GetLatLonAltDouble() { return GetLatLonAltReal().cast<double>(); }
    Real GetLatitudeReal() { return latitude_; }
    Real GetLongitudeReal() { return longitude_; }
    Real GetAltitudeReal() { return altitude_; }
    VecX GetLatLonAltReal() {
      VecX lla(3);
      lla << latitude_, longitude_, altitude_;
      return lla;
    }

    CartesianOrbitState GetCartesianGCRFStateAtEpoch(Real epoch) override {
      // Convert RV to GCRF
      Vec6 rv = GetRvState()->GetVec();
      Vec6 rv_gcrf = ConvertFrame(epoch, rv, frame_, Frame::GCRF);
      return CartesianOrbitState(rv_gcrf, Frame::GCRF);
    }
  };

}  // namespace lupnt
