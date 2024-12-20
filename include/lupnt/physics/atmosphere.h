/**
 * @file atmosphere.h
 * @author Stanford NAV Lab
 * @brief  Atmosphere and ionosphere models
 * @version 0.1
 * @date 2024-12-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <lupnt/core/constants.h>

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace lupnt {

  /** ionosphere */
  double Klobucher(double t_gps, double elevation, double azimuth, double latitude_u,
                   double longitude_u, double freq_Hz, const Vec4d& alpha, const Vec4d& beta);

}  // namespace lupnt
