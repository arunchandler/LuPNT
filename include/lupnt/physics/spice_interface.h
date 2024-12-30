/**
 * @file spice_interface.h
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string.h>

#include <map>
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  namespace spice {
    static segment_t* cheby_s;
    static long cheby_n;
    void LoadSpiceKernel(void);
    void ExtractPckCoeffs(void);
    Mat6d GetFrameConversionMat(Real t_tai, const std::string& from_frame,
                                const std::string& to_frame);
    Vec3d GetPlanetOrientation(NaifId id, Real t_tdb);

    Real String2TDB(const std::string& str);
    Real String2TAI(const std::string& str);

    std::string TAItoStringUTC(Real t_tai, int prec);
    std::string TDBtoStringUTC(Real t_tdb, int prec);

    Real ConvertTime(Real t, Time from_time, Time to_time);

    Vec6 GetBodyPosVel(const Real t_tai, NaifId center, NaifId target);
    MatX6 GetBodyPosVel(const VecX& t_tai, NaifId center, NaifId target);

    Vec3d GetBodyPosSpice(Real t_tai, NaifId obs, NaifId target,
                          const std::string& refFrame = "J2000",
                          const std::string& abCorrection = "NONE");
    Vec6d GetBodyPosVelSpice(Real t_tai, NaifId obs, NaifId target,
                             const std::string& refFrame = "J2000",
                             const std::string& abCorrection = "NONE");
  }  // namespace spice
}  // namespace lupnt
