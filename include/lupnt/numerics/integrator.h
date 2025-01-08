/**
 * @file integrator.h
 * @author Stanford NAV LAB
 * @brief Integrator interfaces
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>

#include "lupnt/core/constants.h"

namespace lupnt {

  using ODE = std::function<VecX(Real, const VecX&)>;

  enum class IntegratorType {
    RK4,
    RK8,
    RKF45,
    PD45,
  };
  constexpr IntegratorType default_integrator = IntegratorType::RK4;

  class IntegratorParams {
  public:
    int max_iter = 20;
    double abstol = 1e-6;
    double reltol = 1e-6;

    IntegratorParams() = default;
    IntegratorParams(int max_iter, double abstol, double reltol)
        : max_iter(max_iter), abstol(abstol), reltol(reltol) {
      CheckIntegratorParams();
    };
    void CheckIntegratorParams();
  };

  class IIntegrator {
  protected:
    IntegratorParams params_;

  public:
    virtual VecX Step(const ODE& f, Real t, const VecX& x, Real dt) = 0;
    void SetIntegratorParams(IntegratorParams params) { params_ = params; };
    virtual ~IIntegrator() {};
  };

  // Runge-Kutta Integrators
  class RK4 : public IIntegrator {
  public:
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt);
  };

  class RK8 : public IIntegrator {
  public:
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt);
  };

  // Runge-Kutta-Fehlberg Integrators with adaptive step size
  class IRKF : public IIntegrator {
  private:
    int order_;

  public:
    IRKF() = default;
    IRKF(IntegratorParams params, int order) : order_(order) { SetIntegratorParams(params); };
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt) override;
    bool ComputeRelError(const VecX& x_new_low, const VecX& x_new_high, Real dt);
    virtual void Update(const ODE& f, Real t, const VecX& x, Real dt, VecX& x_new_low,
                        VecX& x_new_high)
        = 0;
    virtual ~IRKF() = default;
  };

  class RKF45 : public IRKF {
  public:
    RKF45(IntegratorParams params) : IRKF(params, 4) {};
    void Update(const ODE& f, Real t, const VecX& x, Real dt, VecX& x_new_low,
                VecX& x_new_high) override;
  };

  class PD45 : public IIntegrator {
  private:
    static constexpr std::array<std::array<double, 6>, 7> A_
        = {{{0, 0, 0, 0, 0, 0},
            {1.0 / 5, 0, 0, 0, 0, 0},
            {3.0 / 40, 9.0 / 40, 0, 0, 0, 0},
            {44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0},
            {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0, 0},
            {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656, 0},
            {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}}};

    static constexpr std::array<double, 7> b_
        = {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0};

    static constexpr std::array<double, 7> b_star_ = {
        5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40};

  public:
    PD45(IntegratorParams params) { SetIntegratorParams(params); };
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt) override;
  };

}  // namespace lupnt
