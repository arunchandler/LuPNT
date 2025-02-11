/**
 * @file DynamicsNumerical.cpp
 * @author Stanford NAV LAB
 * @brief List of Numerical Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/core/constants.h"
#include "lupnt/core/progress_bar.h"
#include "lupnt/dynamics/dynamics.h"

namespace lupnt {

  NumericalOrbitDynamics::NumericalOrbitDynamics(ODE odefunc, IntegratorType integrator)
      : odefunc_(odefunc), propagator_(integrator) {}

  void NumericalOrbitDynamics::SetTimeStep(Real dt) { dt_ = dt; };
  Real NumericalOrbitDynamics::GetTimeStep() const { return dt_; };
  void NumericalOrbitDynamics::SetODEFunction(ODE odefunc) { odefunc_ = odefunc; };

  Vec6 NumericalOrbitDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    if (abs(tf - t0) < EPS) return x0;
    if (stm == nullptr) {
      Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_);
      return xf;
    } else {
      MatXd stm_X(6, 6);
      Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_, &stm_X);
      *stm = stm_X;
      return xf;
    }
  }

  MatX6 NumericalOrbitDynamics::Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress) {
    MatX6 xf = MatX6::Zero(tf.size(), 6);
    ProgressBar pbar(tf.size());
    pbar.SetDescription("Propagating");
    xf.row(0) = Propagate(x0, t0, tf(0));
    for (int i = 1; i < tf.size(); i++) {
      Real t0_i = tf(i - 1);
      Real tf_i = tf(i);
      Vec6 x0_i = xf.row(i - 1);
      Vec6 xf_i = Propagate(x0_i, t0_i, tf_i);
      xf.row(i) = xf_i;
      if (progress) pbar.Update(i);
    }
    if (progress) pbar.Finish();
    return xf;
  }

  NumericalClockDynamics::NumericalClockDynamics(ODE odefunc, IntegratorType integrator)
      : odefunc_(odefunc), propagator_(integrator) {}

  void NumericalClockDynamics::SetTimeStep(Real dt) { dt_ = dt; }
  Real NumericalClockDynamics::GetTimeStep() const { return dt_; }
  void NumericalClockDynamics::SetPlanetaryStates(MatX planetary_states) { planetary_states_ = planetary_states; }
  MatX NumericalClockDynamics::GetPlanetaryStates() const { return planetary_states_; }
  void NumericalClockDynamics::SetODEFunction(ODE odefunc) { odefunc_ = odefunc; }

  VecX NumericalClockDynamics::Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm) {
    if (abs(tf - t0) < EPS) return x0;
    if (stm == nullptr) {
      VecX xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_);
      return xf;
    } else {
      MatXd stm_X(6, 6);
      VecX xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_, &stm_X);
      *stm = stm_X;
      return xf;
    }
  }

  // ****************************************************************************
  // CircularRestrictedThreeBodyDynamics
  // ****************************************************************************

  CR3BPDynamics::CR3BPDynamics(Real mu, IntegratorType integ)
    : NumericalOrbitDynamics([this](Real t, const Vec6 &x) { return ComputeRates(t, x); }, integ),
      mu_(mu) {}

  Vec6 CR3BPDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Real x_val = x(0);
    Real y_val = x(1);
    Real z_val = x(2);
    Real x_dot = x(3);
    Real y_dot = x(4);
    Real z_dot = x(5);

    Real r1 = sqrt(pow(x_val + mu_, 2) + pow(y_val, 2) + pow(z_val, 2));
    Real r2 = sqrt(pow(x_val - 1 + mu_, 2) + pow(y_val, 2) + pow(z_val, 2));

    Real ax = 2 * y_dot + x_val - (1 - mu_) * (x_val + mu_) / pow(r1, 3) - mu_ * (x_val - 1 + mu_) / pow(r2, 3);
    Real ay = -2 * x_dot + y_val - (1 - mu_) * y_val / pow(r1, 3) - mu_ * y_val / pow(r2, 3);
    Real az = -(1 - mu_) * z_val / pow(r1, 3) - mu_ * z_val / pow(r2, 3);

    Vec6 rates;
    rates << x_dot, y_dot, z_dot, ax, ay, az;

    return rates;
  }

  OrbitState CR3BPDynamics::PropagateState(const OrbitState &state, Real t0, Real tf, Mat6d *stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CARTESIAN);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  // ****************************************************************************
  // CartesianTwoBodyDynamics
  // ****************************************************************************

  CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(Real GM, IntegratorType integ)
      : NumericalOrbitDynamics([this](Real t, const Vec6 &x) { return ComputeRates(t, x); }, integ),
        GM_(GM) {};

  Vec6 CartesianTwoBodyDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Vec6 rv_dot;
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    return rv_dot;
  }

  OrbitState CartesianTwoBodyDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                      Mat6d *stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CARTESIAN);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  // ****************************************************************************
  // J2CartTwoBodyDynamics
  // ****************************************************************************

  J2CartTwoBodyDynamics::J2CartTwoBodyDynamics(Real GM, Real J2, Real R_body, IntegratorType integ)
      : NumericalOrbitDynamics([this](Real t, const Vec6 &x) { return ComputeRates(t, x); }, integ),
        GM_(GM),
        J2_(J2),
        R_body_(R_body) {};

  Vec6 J2CartTwoBodyDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Vec6 rv_dot(6);
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    Vec3 a_J2;
    Real aux1 = -3.0 / 2.0 * GM_ * J2_ * pow(R_body_, 2.0) / pow(r_norm, 4.0);
    Real aux2 = 5.0 * pow(r(2) / R_body_, 2.0);
    a_J2(0) = aux1 * (1.0 - aux2) * r(0);
    a_J2(1) = aux1 * (1.0 - aux2) * r(1);
    a_J2(2) = aux1 * (3.0 - aux2) * r(2);

    rv_dot.tail(3) += a_J2;

    return rv_dot;
  }

  OrbitState J2CartTwoBodyDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                   Mat6d *stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CARTESIAN);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  // ****************************************************************************
  // J2KeplerianDynamics
  // ****************************************************************************

  J2KeplerianDynamics::J2KeplerianDynamics(Real GM, Real J2, Real R_body, IntegratorType integ)
      : NumericalOrbitDynamics([this](Real t, const Vec6 &x) { return ComputeRates(t, x); }, integ),
        GM_(GM),
        J2_(J2),
        R_body_(R_body) {};

  Vec6 J2KeplerianDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Real p = x(0) * (1.0 - x(1) * x(1));
    Real n = sqrt(GM_ / pow(x(0), 3.0));
    Real eta = sqrt(1.0 - x(1) * x(1));

    Vec6 coe_dot(6);
    coe_dot(0) = 0.0;
    coe_dot(1) = 0.0;
    coe_dot(2) = 0.0;
    coe_dot(3) = -3.0 / 2.0 * J2_ * pow(R_body_ / p, 2.0) * n * cos(x(2));
    coe_dot(4) = 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * (5.0 * pow(cos(x(2)), 2.0) - 1.0);
    coe_dot(5)
        = n + 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * eta * (3.0 * pow(cos(x(2)), 2.0) - 1.0);
    return coe_dot;
  }

  OrbitState J2KeplerianDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                 Mat6d *stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CLASSICAL_OE);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  // ****************************************************************************
  // MoonMeanDynamics
  // ****************************************************************************

  MoonMeanDynamics::MoonMeanDynamics(IntegratorType integ)
      : NumericalOrbitDynamics([this](Real t, const Vec6 &x) { return ComputeRates(t, x); },
                               integ) {};

  Vec6 MoonMeanDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Real a = x(0);
    Real e = x(1);
    Real i = x(2);
    Real w = x(4);

    Vec6 coe_dot(6);
    coe_dot(0) = 0;
    coe_dot(1) = (15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON)) * e
                 * sqrt(1 - pow(e, 2)) * pow(sin(i), 2) * sin(2 * w);
    coe_dot(2) = -(15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (16 * sqrt(GM_MOON)) * pow(e, 2)
                 / sqrt(1 - pow(e, 2)) * sin(2 * i) * sin(2 * w);
    coe_dot(3) = -(3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (2 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * cos(i)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * (5 * pow(e, 2) * cos(2 * w) - 3 * pow(e, 2) - 2) * cos(i);
    coe_dot(4) = (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * (5 * pow(cos(i), 2) - 1)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * ((5 * pow(cos(i), 2) - 1 + pow(e, 2))
                          + 5 * (1 - pow(e, 2) - pow(cos(i), 2)) * cos(2 * w));
    coe_dot(5) = sqrt(GM_MOON) / pow(a, 3.0 / 2.0)
                 + (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                       / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 3.0 / 2.0))
                       * (3 * pow(cos(i), 2) - 1)
                 - (k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON))
                       * ((3 * pow(e, 2) + 7) * (3 * pow(cos(i), 2) - 1)
                          + 15 * (1 + pow(e, 2)) * pow(sin(i), 2) * cos(2 * w));
    return coe_dot;
  }

  OrbitState MoonMeanDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                              Mat6d *stm) {
    CheckOrbitStateRepres(state, OrbitStateRepres::CLASSICAL_OE);
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

// ****************************************************************************
// RelativityClockDynamics
// ****************************************************************************

RelativityClockDynamics::RelativityClockDynamics(IntegratorType integ)
    : NumericalClockDynamics([this](Real t, const VecX &x) { return ComputeRates(t, x); }, integ) {}

VecX RelativityClockDynamics::ComputeRates(Real t, const VecX &x) const {

    Vec3 pos = x.segment(0, 3);
    Vec3 vel = x.segment(3, 3);

    MatX planetary_states = GetPlanetaryStates();
    Real U = 0;

    // Compute the gravitational potential difference U
    Real r = pos.norm();
    for (int i = 0; i < planetary_states.cols(); ++i) {
      Vec3 planet_pos = planetary_states.block<3, 1>(0, i);
      Real planet_mu = planetary_states(6, i);
      Real r_planet = (pos - planet_pos).norm();
      U += -planet_mu / r_planet;
    }
    // Compute squared velocity
    Real V2 = vel.squaredNorm();

    // Compute relativistic clock rate correction
    Real dT_dtau = 1 + (U / (C * C)) + (0.5 * V2 / (C * C));

    VecX rates;
    rates << dT_dtau;

    return rates;
  }

};  // namespace lupnt
