/**
 * @file clock.h
 * @author Stanford NAV LAB
 * @brief Clock class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <tuple>

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/numerics/math_utils.h"

#include "state.h"

namespace lupnt {

  enum class ClockModel { kMicrosemiCsac, kRafs, kUso, kMiniRafs, kUnknown };

  std::tuple<Real, Real, Real> GetClockSigma(ClockModel clk_model);

  class ClockState : public IState {
  private:
    VecX x_;
    int state_size_;  // it can be 2 or 3

  public:
    ClockState();
    ClockState(int state_size);
    ClockState(VecX clock_vec);

    // Overrides
    int GetSize() const override;
    VecX GetVec() const override;
    void SetVec(const VecX& x) override;
    Real GetValue(int i) const override;
    void SetValue(int idx, Real val) override;
    StateType GetStateType() const override;
  };

  class ClockDynamics : public IDynamics {
  private:
    ClockModel clk_model_ = ClockModel::kUnknown;
    bool noise_ = false;

  public:
    ClockDynamics();
    ClockDynamics(ClockModel clk_model);

    inline void SetNoise(bool noise) { noise_ = noise; }

    static Mat2 TwoStatePhi(Real dt);
    static Mat3 ThreeStatePhi(Real dt);

    static Mat2 TwoStateNoise(ClockModel clk_model, Real dt);
    static Mat3 ThreeStateNoise(ClockModel clk_model, Real dt);

    using IDynamics::Propagate;
    ClockState PropagateState(const ClockState state, Real t0, Real tf, MatXd* stm = nullptr);
    Ptr<IState> PropagateState(const Ptr<IState>& state, Real t0, Real tf,
                               MatXd* stm = nullptr) override;
    VecX Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm = nullptr) override;
    Vec2 Propagate(const Vec2& x0, Real t0, Real tf, Mat2* stm = nullptr);
    Vec3 Propagate(const Vec3& x0, Real t0, Real tf, Mat3* stm = nullptr);
  };

  // Numerical Clock Dynamics Interface
  class NumericalClockDynamics : public ClockDynamics {
    private:
      ODE odefunc_;
      NumericalPropagator propagator_;
      Real dt_ = 10.0;
  
    public:
      NumericalClockDynamics(ODE odefunc = nullptr, IntegratorType integ = default_integrator);
      void SetTimeStep(Real dt);
      Real GetTimeStep() const;
      void SetODEFunction(ODE odefunc);
      void SetIntegratorParams(IntegratorParams params) {
        propagator_.integrator->SetIntegratorParams(params);
      }
  
      // Overrides
      VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) override;
  
      // Interface
      virtual VecX ComputeRates(Real t, const VecX &x) const = 0;
  };

  // Relativistic Clock Dynamics Interface
  template <typename T = double> class RelativisticClockDynamics : public NumericalClockDynamics {
    private:
      std::vector<BodyT<T>> bodies_;
      Frame frame_ = Frame::NONE;
      Vec6 observer_state_;

    public:
      RelativisticClockDynamics(IntegratorType integ = default_integrator);

      void AddBody(const BodyT<T> &body) {
        for (auto &b : bodies_) {
          if (b.id == body.id) throw std::runtime_error("Body already added");
        }
        bodies_.push_back(body);
      }
  
      std::vector<BodyT<T>> GetBodies() { return bodies_; }
  
      void RemoveBody(const BodyT<T> &body) {
        for (auto it = bodies_.begin(); it != bodies_.end(); ++it) {
          if (it->id == body.id) {
            bodies_.erase(it);
            break;
          }
        }
      }

      void SetFrame(Frame frame) { frame_ = frame; }
      void GetFrame(Frame &frame) { frame = frame_; }

      void SetObserverState(const Vec6 &state) { observer_state_ = state; }
      Vec6 GetObserverState() { return observer_state_; }

      // Overrides
      VecX ComputeRates(Real t, const VecX &x) const override;
  };

  class ClockOrbitDynamics : public IDynamics {
    private:
      ODE odefunc_;
      NumericalPropagator propagator_;
      Real dt_ = 10.0;

    public:
      std::shared_ptr<NBodyDynamics<>> orbitDynamics_;
      std::shared_ptr<RelativisticClockDynamics<>> clockDynamics_;
      ClockOrbitDynamics(IntegratorType integ = default_integrator);

      void AddBody_COD(const BodyT<> &body) {
        orbitDynamics_->AddBody(body);
        clockDynamics_->AddBody(body);
      }

      void RemoveBody_COD(const BodyT<> &body) {
        orbitDynamics_->RemoveBody(body);
        clockDynamics_->RemoveBody(body);
      }

      void SetTimeStep(Real dt) { dt_ = dt; orbitDynamics_->SetTimeStep(dt); clockDynamics_->SetTimeStep(dt); }

      void SetFrame_COD(const Frame &frame) { orbitDynamics_->SetFrame(frame); clockDynamics_->SetFrame(frame); }

      std::vector<BodyT<>> GetBodies_COD() {
        std::vector<BodyT<>> bodies = orbitDynamics_->GetBodies();
        return bodies;
      }

      void SetObserverState_COD(const Vec6 &state) { clockDynamics_->SetObserverState(state); }

      VecX ComputeRates(Real t, const VecX &x) const;

      VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr);

      //TODO: Implement this
      Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf, MatXd *stm = nullptr);

  };

}  // namespace lupnt
