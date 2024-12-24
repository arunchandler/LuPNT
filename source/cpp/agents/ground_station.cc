/**
 * @file groundstation.cc
 * @author Stanford NAV Lab
 * @brief  GroundStation Agent
 * @version 0.1
 * @date 2024-12-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <lupnt/agents/ground_station.h>
#include <lupnt/core/file.h>

#include <fstream>

namespace lupnt {

  GroundStationDataMap LoadGroundStationData() {
    // Path to the ground station data file
    std::string gs_data_file = GetGroundStationDataDir() / "ground_stations.csv";
    std::ifstream file = OpenFile<std::ifstream>(gs_data_file);
    std::string line;
    std::getline(file, line);  // Skip the header
    // gs_id, name, location, dish_size, lat_deg, lon_deg, alt_m, eirp_min, eirp_max,
    // uplink_freq_min, uplink_freq_max, downlink_freq_min, downlink_freq_max,
    // gt_dBk, A, D, R, freq_band

    int gs_idx = 0;
    GroundStationDataMap ground_stations_data;

    while (std::getline(file, line)) {
      // skip the header line
      if (line.find("gs_id") != std::string::npos) continue;

      std::stringstream ss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (std::getline(ss, token, ',')) {
        tokens.push_back(token);
      }

      GroundStationData gs_data;
      gs_data.gs_id = gs_idx;
      gs_data.name = tokens[1];
      gs_data.location = tokens[2];
      gs_data.dish_size = std::stod(tokens[3]);
      gs_data.latitude = std::stod(tokens[4]);
      gs_data.longitude = std::stod(tokens[5]);
      gs_data.altitude_m = std::stod(tokens[6]);
      gs_data.eirp_min = std::stod(tokens[7]);
      gs_data.eirp_max = std::stod(tokens[8]);
      gs_data.uplink_freq_min = std::stod(tokens[9]);
      gs_data.uplink_freq_max = std::stod(tokens[10]);
      gs_data.downlink_freq_min = std::stod(tokens[11]);
      gs_data.downlink_freq_max = std::stod(tokens[12]);
      gs_data.gt_dBk = std::stod(tokens[13]);
      gs_data.angle_meas = std::stoi(tokens[14]);
      gs_data.doppler_meas = std::stoi(tokens[15]);
      gs_data.range_meas = std::stoi(tokens[16]);
      gs_data.freq_band = tokens[17];

      ground_stations_data[gs_data.gs_id] = gs_data;

      gs_idx += 1;
    }

    std::cout << "Loaded " << gs_idx << " ground stations" << std::endl;

    return ground_stations_data;
  }

}  // namespace lupnt
