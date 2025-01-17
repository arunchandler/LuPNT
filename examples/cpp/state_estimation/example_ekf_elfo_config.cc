/**
 * @file example_ekf_elfo_config.cc
 * @author Stanford NAVLAB
 * @brief lunar orbit estimation using GPS measurements with yaml config
 * @version 0.1
 * @date 2024-12-28
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  // Read config file
  std::string config_file = "ex_nav_elfo_gps.yaml";
  ConfigReader reader(config_file, true);

  // Initialize simulation
  // NavSimulation sim;

  return 0;
}