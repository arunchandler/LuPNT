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
    Real t_diff_in = ConvertTime(t_tai, Time::TAI, time_in);
    for (auto time_out : time_list) {
      Real t_out = ConvertTime(t_diff_in, time_in, time_out);
      Real t_final = ConvertTime(t_out, time_out, Time::TAI);
      RequireNear(t_tai, t_final, epsilon);
    }
  }

  // Spice
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

  // Orekit
  auto file = OpenTestDataFile("time_converter_orekit.txt");
  std::string line;
  int precision = 5;
  Real t_utc;
  while (std::getline(file, line)) {
    std::string time, date_in;
    double mjd_in, t_diff_in;
    std::istringstream iss(line);
    iss >> time >> date_in >> mjd_in >> t_diff_in;
    if (time == "UTC") t_utc = Gregorian2Time(date_in);
    Real t = ConvertTime(t_utc, Time::UTC, string2time.at(time));
    Real t_j2000 = ConvertTime(0, Time::TT, string2time.at(time));
    Real t_diff = t - t_j2000;
    Real mjd = Time2MJD(t);
    std::string date = Time2GregorianString(t, precision);
    RequireNear(mjd, mjd_in, epsilon);
    RequireNear(t_diff, t_diff_in, epsilon);
    int len = date_in.size();
    if (date.size() < len) len = date.size();
    len--;
    REQUIRE_THAT(date.substr(0, len), Equals(date_in.substr(0, len)));
  }
}
