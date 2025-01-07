/**
 * @file FrameConverter.cpp
 * @author Stanford NAV LAB
 * @brief Coordinate conversion functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <ostream>

#include "cheby.h"
#include "lupnt/core/constants.h"

namespace lupnt {

  class CartesianOrbitState;

  enum class Frame {
    // Earth
    NONE,         // No frame
    ICRF,         // International Celestial Reference System
    ITRF,         // International Terrestrial Reference Frame
    ECEF = ITRF,  // Earth-Centered Earth-Fixed
    GCRF,         // Geocentric Reference System
    EME,          // Earth-Centered mean equator and equinox at J2000
    ECI = EME,    // Earth-Centered Inertial
    SER,          // Sun-Earth Rotating Frame
    GSE,          // Geocentric Solar Ecliptic
    MOD,          // Mean of date equatorial system
    TOD,          // True of date equatorial system
    EMR,          // Earth-Moon Rotating Frame
    // Moon
    MOON_CI,  // Moon-centered Inertial Frame (Axis aligened with ICRF)
    MOON_PA,  // Moon-Fixed with principal axes
    MOON_ME,  // Moon-Fixed with mean-Earth / polar axes
    MOON_OP,  // Earth Orbit Frame
    // Solar System
    MERCURY_FIXED,  // Mercury fixed frame
    VENUS_FIXED,    // Venus fixed frame
    MARS_FIXED,     // Mars fixed frame
    JUPITER_FIXED,  // Jupiter fixed frame
    SATURN_FIXED,   // Saturn fixed frame
    URANUS_FIXED,   // Uranus fixed frame
    NEPTUNE_FIXED,  // Neptune fixed frame
    // Inertial
    MERCURY_CI,  // Mercury-centered Inertial Frame
    VENUS_CI,    // Venus-centered Inertial Frame
    MARS_CI,     // Mars-centered Inertial Frame
    JUPITER_CI,  // Jupiter-centered Inertial Frame
    SATURN_CI,   // Saturn-centered Inertial Frame
    URANUS_CI,   // Uranus-centered Inertial Frame
    NEPTUNE_CI,  // Neptune-centered Inertial Frame
  };

  std::ostream &operator<<(std::ostream &os, Frame frame);

  extern const std::map<Frame, NaifId> frame_centers;

  extern std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6 &rv)>>
      frame_conversions;

  // Vec = func(real, Vec)
  Vec6 ConvertFrame(Real t_tai, const Vec6 &rv_in, Frame frame_in, Frame frame_out,
                    bool rotate_only = false);
  Vec3 ConvertFrame(Real t_tai, const Vec3 &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(real, Mat)
  MatX6 ConvertFrame(Real t_tai, const MatX6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);
  MatX3 ConvertFrame(Real t_tai, const MatX3 &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(Vec, Vec)
  MatX6 ConvertFrame(VecX t_tai, const Vec6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);
  MatX3 ConvertFrame(VecX t_tai, const Vec3 &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(Vec, Mat)
  MatX6 ConvertFrame(VecX t_tai, const MatX6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);
  MatX3 ConvertFrame(VecX t_tai, const MatX3 &r_in, Frame frame_in, Frame frame_out);

  CartesianOrbitState ConvertFrame(Real t_tai, const CartesianOrbitState &state_in, Frame frame_out,
                                   bool rotate_only = false);

}  // namespace lupnt
