#include <lupnt/data/kernels.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

// Test cases for time conversion
TEST_CASE("data.Kernels") {
  // Spice
  file = OpenTestDataFile("body_pos_vel_cspice.txt");
  while (std::getline(file, line)) {
    std::string time, date_in;
    double t_diff_in;
    std::istringstream iss(line);
    iss >> time >> date_in >> t_diff_in;
    if (time == "UTC") t_utc = Gregorian2Time(date_in);
    Real t = ConvertTime(t_utc, Time::UTC, string2time.at(time));
    Real t_j2000 = ConvertTime(0, Time::TT, string2time.at(time));
    Real t_diff = t - t_j2000;
    RequireNear(t_diff, t_diff_in, epsilon);
  }
}
