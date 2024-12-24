#include <lupnt/physics/spice_interface.h>
#include <lupnt/physics/time_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

std::vector<Time> time_list
    = {Time::UT1, Time::UTC, Time::TAI, Time::TDB, Time::TT, Time::TCG, Time::GPS};

// Test cases for time conversion
TEST_CASE("physics.ConvertTime") {
  // All-to-all time conversions
  for (auto time_in : time_list) {
    Real t_tai = Gregorian2Time(2022, 12, 18, 18, 0, 0.0);
    Real t_in = ConvertTime(t_tai, Time::TAI, time_in);
    for (auto time_out : time_list) {
      Real t_out = ConvertTime(t_in, time_in, time_out);
      Real t_final = ConvertTime(t_out, time_out, Time::TAI);
      RequireNear(t_tai, t_final, epsilon);
    }
  }

  std::vector<std::vector<double>> dates = {{2000, 1, 1, 12, 0, 0.0}, {2022, 12, 18, 18, 0, 0.0}};
  std::vector<std::string> dates_str = {"2000 Jan 01, 12:00:00 TDB", "2022 Dec 18, 18:00:00 TDB"};
  for (int i = 0; i < dates.size(); i++) {
    Real t_tdb = Gregorian2Time(dates[i][0], dates[i][1], dates[i][2], dates[i][3], dates[i][4],
                                dates[i][5]);
    Real t_tdb_sp = spice::String2TDB(dates_str[i]);
    RequireNear(t_tdb, t_tdb_sp, epsilon);

    Real t_tai = ConvertTime(t_tdb, Time::TDB, Time::TAI);
    Real t_tai_sp = spice::ConvertTime(t_tdb, Time::TDB, Time::TAI);
    RequireNear(t_tai, t_tai_sp, epsilon);

    Real t_tt = ConvertTime(t_tai, Time::TAI, Time::TT);
    Real t_tt_sp = spice::ConvertTime(t_tai, Time::TAI, Time::TT);
    RequireNear(t_tt, t_tt_sp, epsilon);
  }
}
