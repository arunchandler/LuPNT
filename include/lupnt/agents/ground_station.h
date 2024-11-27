/**
 * @file ground_stations.h
 * @author your name (you@domain.com)
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

enum GsID {
    NONE,
    DSS13,
    DSS14,
    DSS15,
    DSS24,
    DSS25,
    DSS26,
    DSS34,
    DSS35,
    DSS36,
    DSS43,
    DSS45,
    DSS54,
    DSS55,
    DSS63,
    DSS65,
    SSC_Clewistion,
    SSC_Esrange,
    SSC_Inuvik,
    SSC_Northpole,
    SSC_PuntaArenas,
    SSC_Santiago,
    SSC_Siracha,
    SSC_SouthPoint,
    SSC_WASC,
    SSC_Stockholm,
    LEGS_WhiteSands,
    LEGS_SouthAfrica,
    ESTRACK_CEB,
    ESTRACK_KI1,
    ESTRACK_KI2,
    ESTRACK_KRU,
    ESTRACK_MLG,
    ESTRACK_NNO,
    ESTRACK_NNO2,
    ESTRACK_RED,
    ESTRACK_SMA,
};

  struct GroundStationData {
    GsID gs_id;
    std::string name;
    double latitude;
    double longitude;
    double altitude_m;
  };

  std::vector<GroundStationData> GROUND_STATION_LIST = {
      // DSN:
      // Reference: https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf
      // Goldstone
      {DSS13, "DSS_13", DegMinSecToRad(35, 14, 49.79131), DegMinSecToRad(243, 12, 19.94761), 1070.444},
      {DSS14, "DSS_14", DegMinSecToRad(35, 25, 33.24312), DegMinSecToRad(243, 6, 37.66244), 1001.390},
      {DSS15, "DSS_15", DegMinSecToRad(35, 25, 18.67179), DegMinSecToRad(243, 6, 46.09762), 973.211},
      {DSS24, "DSS_24", DegMinSecToRad(35, 20, 23.61416), DegMinSecToRad(243, 7, 30.74007), 951.499},
      {DSS25, "DSS_25", DegMinSecToRad(35, 20, 15.40306), DegMinSecToRad(243, 7, 28.69246), 959.634},
      {DSS26, "DSS_26", DegMinSecToRad(35, 20, 8.48118), DegMinSecToRad(243, 7, 37.14062), 968.686},
      // Camberra
      {DSS34, "DSS_34", DegMinSecToRad(-35, 23, 54.52383), DegMinSecToRad(148, 58, 55.07191), 692.020},
      {DSS35, "DSS_35", DegMinSecToRad(-35, 23, 44.86387), DegMinSecToRad(148, 58, 53.24088), 694.897},
      {DSS36, "DSS_36", DegMinSecToRad(-35, 23, 42.36634), DegMinSecToRad(148, 58, 42.75912), 685.503},
      {DSS43, "DSS_43", DegMinSecToRad(-35, 24, 8.72724), DegMinSecToRad(148, 58, 52.56231), 688.867},
      {DSS45, "DSS_45", DegMinSecToRad(-35, 23, 54.44766), DegMinSecToRad(148, 58, 39.66828), 674.347},
      // Madrid
      {DSS54, "DSS_54", DegMinSecToRad(40, 25, 32.23805), DegMinSecToRad(355, 44, 45.25141), 837.051},
      {DSS55, "DSS_55", DegMinSecToRad(40, 25, 27.46525), DegMinSecToRad(355, 44, 50.52012), 819.061},
      {DSS63, "DSS_63", DegMinSecToRad(40, 25, 52.35510), DegMinSecToRad(355, 45, 7.16924), 864.816},
      {DSS65, "DSS_65", DegMinSecToRad(40, 25, 37.94289), DegMinSecToRad(355, 44, 57.48397), 833.854},
      // SSC stations
      // Reference: https://sscspace.com/services/satellite-ground-stations/our-stations/
      // Todo: Convert altitude into WGS84
      {SSC_Clewistion, "SSC_Clewistion", DegMinSecToRad(26, 44, 0.0), DegMinSecToRad(-81, 2, 0.0), 6.0},
      {SSC_Esrange, "SSC_Esrange", DegMinSecToRad(67, 53, 0.0), DegMinSecToRad(21, 4, 0.0), 390.0},
      {SSC_Inuvik, "SSC_Inuvik", DegMinSecToRad(68, 24, 0.0), DegMinSecToRad(-133, 30, 0.0), 88.0},
      {SSC_Northpole, "SSC_Northpole", DegMinSecToRad(64, 48, 0.0), DegMinSecToRad(-147, 30, 0.0), 142.0},
      {SSC_PuntaArenas, "SSC_PuntaArenas", DegMinSecToRad(-52, 56, 0.0), DegMinSecToRad(-70, 51, 0.0), 19.0},
      {SSC_Santiago, "SSC_Santiago", DegMinSecToRad(-33, 8, 0.0), DegMinSecToRad(-70, 40, 0.0), 694.0},
      {SSC_Siracha, "SSC_Siracha", DegMinSecToRad(13, 6, 0.0), DegMinSecToRad(100, 55, 0.0), 20.0},
      {SSC_SouthPoint, "SSC_SouthPoint", DegMinSecToRad(19, 1, 0.0), DegMinSecToRad(-155, 40, 0.0), 356.0},
      {SSC_WASC, "SSC_WASC", DegMinSecToRad(-29, 5, 0.0), DegMinSecToRad(115, 35, 0.0), 281.0},
      {SSC_Stockholm, "SSC_Stockholm", DegMinSecToRad(59, 21, 0.0), DegMinSecToRad(18, 5, 0.0), 15.00},
      // LEGS
      // https://explorers.larc.nasa.gov/2023ESE/pdf_files/LEGS%20Brochure%20r20.pdf
      // The third legs point is not available yet
      // Todo: Convert altitude into WGS84
      {LEGS_WhiteSands, "LEGS_WhiteSands", DegMinSecToRad(32, 32, 41.5), DegMinSecToRad(-106, 36, 45.0), 1462.0},
      {LEGS_SouthAfrica, "LEGS_SouthAfrica", DegMinSecToRad(-33, 13, 52.4), DegMinSecToRad(20, 34, 53.9), 900.0},
      // ESA Stations
      // http://estracknow.esa.int/#/2024-10
      {ESTRACK_CEB, "ESTRACK_CEB", 40.4526901 * RAD, -4.3675499 * RAD, 794.09},
      {ESTRACK_KI1, "ESTRACK_KI1", 67.8571243 * RAD, 20.964325 * RAD, 402.17},
      {ESTRACK_KI2, "ESTRACK_KI2", 67.858429 * RAD, 20.9668808 * RAD, 400.68},
      {ESTRACK_KRU, "ESTRACK_KRU", 5.2514391 * RAD, -52.8046646 * RAD, -14.67},
      {ESTRACK_MLG, "ESTRACK_MLG", -35.7760086 * RAD, -69.3981934 * RAD, 1550.0},
      {ESTRACK_NNO, "ESTRACK_NNO", -31.0482254 * RAD, 116.1914978 * RAD, 252.26},
      {ESTRACK_NNO2, "ESTRACK_NNO2", -31.0488949 * RAD, 116.1888123 * RAD, 262.94},
      {ESTRACK_RED, "ESTRACK_RED", 50.0004578* RAD, 5.1453438 * RAD, 368.68},
      {ESTRACK_SMA, "ESTRACK_SMA", 36.9972496* RAD, -25.1357212 * RAD, 275.00},

  };

  GroundStationData GetGroundStationData(GsID gs_id) {
    for (auto gs : GROUND_STATION_LIST) {
      if (gs.gs_id == gs_id) {
        return gs;
      }
    }
    return GroundStationData();
  }


  class GroundStation : public Agent {

    protected:
      GsID id_ = GsID::NONE;
      std::string name_ = "NONE";
      NaifId body_id_ = NaifId::EARTH;
      Real latitude_ = 0.0;   // [rad]
      Real longitude_ = 0.0;  // [rad]
      Real altitude_ = 0.0;   // [km]

    public:
      GroundStation(GroundStationData gs_data, NaifId body_id) : Agent() { 
          id_ = gs_data.gs_id;
          name_ = gs_data.name;
          latitude_ = gs_data.latitude;
          longitude_ = gs_data.longitude;
          altitude_ = gs_data.altitude_m * 1e-3;

          SetIsBodyFixed(true); 
          SetBodyId(body_id_);
          // Create body for corresponding body_id
          BodyData body_data = GetBodyData(body_id_);

          SetDynamics(std::make_shared<SurfaceStaticDynamics>(body_id_, body_data.fixed_frame));

          // convert (lat, lon) to ECEF
          Vec3 lla = {latitude_, longitude_, altitude_};
          Vec3 pos_f = LatLonAlt2Cart(lla, body_data.R, body_data.flattening);

          // Set Cartesian State
          Frame frame = body_data.fixed_frame;
          Vec6 rv = {pos_f(0), pos_f(1), pos_f(2), 0.0, 0.0, 0.0};
          SetRvState(std::make_shared<CartesianOrbitState>(rv, frame));
      };
    };
} // namespace lupnt
