/**
 * @file NBodyDynamics.cpp
 * @author Stanford NAV LAB
 * @brief Multiple-body dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/core/constants.h"
#include "lupnt/data/kernels.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/dynamics/forces.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/solar_system.h"

namespace lupnt {

  template <typename T> NBodyDynamics<T>::NBodyDynamics(IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&NBodyDynamics::ComputeRates, this, std::placeholders::_1,
                                         std::placeholders::_2),
                               integ){};
  template class NBodyDynamics<double>;
  template class NBodyDynamics<Real>;

  template <typename T> Vec6 NBodyDynamics<T>::ComputeRates(Real t_tai, const Vec6& rv) const {
    if (frame_ == Frame::NONE) throw std::runtime_error("Frame not set");

    // Position, velocity, and acceleration [km, km/s, km/s^2]
    // w.r.t. to the inertial frame origin
    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);
    Vec3 a = Vec3::Zero();

    bool only_earth_moon = true; //Earth and Moon CI have the same frame + offset
    // Checking if only Earth and Moon are used
    for (const auto& body : bodies_) {
      if (body.id != NaifId::EARTH && body.id != NaifId::MOON) {
          only_earth_moon = false;
          break;
      }
    }

    for (const auto& body : bodies_) {
      if (body.use_gravity_field) {

        auto& grav = body.gravity_field;
        Vec3 ai;

        if (only_earth_moon) {
          // Position (body-fixed) [km]
          Vector<T, 3> r_bf = ConvertFrame(t_tai, r, frame_, body.fixed_frame).template cast<T>();
          // Acceleration (body-fixed) [km/s^2]
          Vec3 a_bf = AccelarationGravityField<T>(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
          // Acceleration (inertial) [km/s^2]
          ai = ConvertFrame(t_tai, a_bf, body.fixed_frame, body.inertial_frame);
        } else {
          Vec3 r_body = GetBodyPos(t_tai, body.id, frame_);
          ai = AccelerationPointMass(rv.head(3), r_body, body.GM);
          //std::cerr << "Warning: Using only the first order term for acceleration computation for bodies outside of Earth and Moon due to framing issue.\n";
        }
        a += ai;
      } else {
        // Body position w.r.t. the inertial frame origin [km]
        Vec3 r_body = GetBodyPos(t_tai, body.id, frame_);
        // Acceleration (inertial) [km/s^2]
        Vec3 ai = AccelerationPointMass(rv.head(3), r_body, body.GM);
        a += ai;
      }

      // Solar radiation pressure
      if (use_srp_ && body.id != NaifId::SUN) {
        Vec3 r_sun = GetBodyPos(t_tai, body.id, NaifId::SUN, frame_);
        Vec3 a_srp = Illumination(r, r_sun, body.R)
                     * AccelerationSolarRadiation(r, r_sun, area_, mass_, CR_, P_SUN, AU);
        a += a_srp;
      }

      // Atmospheric drag
      if (use_drag_ && body.id == NaifId::EARTH) {
        // TODO: Currently only works for Earth
        Real tt = ConvertTime(t_tai, Time::TAI, Time::TT);
        Real mjd_tt = Time2MJD(tt);
        MatX3 Rot = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000_TT, mjd_tt);
        Vec3 a_drag = AccelerationDrag(mjd_tt, rv, Rot, area_, mass_, CD_);
        a += a_drag;
      }
    }

    Vec6 rv_dot;
    rv_dot << v, a;
    return rv_dot;
  }

  template <typename T> OrbitState NBodyDynamics<T>::PropagateState(const OrbitState& state,
                                                                    Real t0, Real tf, Mat6d* stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CARTESIAN);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

}  // namespace lupnt
