#include <lupnt/core/constants.h>
#include <lupnt/data/kernels.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <magic_enum/magic_enum.hpp>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

// Test cases for time conversion
TEST_CASE("data.Kernels") {
  // Spice
  auto file = OpenTestDataFile("body_pos_vel_cspice.txt");
  std::string line;

  // Time
  std::getline(file, line);
  std::istringstream iss(line);
  std::string tmp;
  double t_tai;
  iss >> tmp >> t_tai;

  // Header
  std::getline(file, line);

  // Values
  while (std::getline(file, line)) {
    std::string center_str, target_str;
    double lt, x, y, z, vx, vy, vz;
    iss.clear();
    iss.str(line);
    iss >> center_str >> target_str >> lt >> x >> y >> z >> vx >> vy >> vz;
    Vec6 rv_in{x, y, z, vx, vy, vz};
    auto center_tmp = enum_cast<NaifId>(center_str);
    auto target_tmp = enum_cast<NaifId>(target_str);
    auto center = enum_cast<NaifId>(center_str).value();
    auto target = enum_cast<NaifId>(target_str).value();
    Vec6 rv = GetBodyPosVel(t_tai, center, target, Frame::GCRF);
    RequireNear(rv, rv_in, epsilon);
  }
}
