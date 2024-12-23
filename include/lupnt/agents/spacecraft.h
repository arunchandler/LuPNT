/**
 * @file spacecraft.h
 * @author Stanford NAV Lab
 * @brief  Spacecraft Agent
 * @version 0.1
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <lupnt/agents/agent.h>

namespace lupnt {
  /**
   * @brief Spacecraft Agent
   *
   */
  class Spacecraft : public Agent {
  protected:
    Ptr<OrbitState> orbit_state_;
    bool orbit_state_set_ = false;

  public:
    Spacecraft() : Agent() { SetIsBodyFixed(false); };

    void SetOrbitState(Ptr<OrbitState> orbit_state) {
      orbit_state_ = orbit_state;
      std::shared_ptr<IState> state = std::static_pointer_cast<OrbitState>(orbit_state);
      Agent::SetRvState(state);
      orbit_state_set_ = true;
    }

    Ptr<OrbitState> GetOrbitState() const {
      // First Update the orbit state vector using state_ vector
      if (!orbit_state_set_) {
        std::cerr << "Orbit State is not set: Call SetOrbitState(Ptr<OrbitState>)" << std::endl;
      }
      orbit_state_->SetVec(rv_->GetVec());
      return orbit_state_;
    }

    CartesianOrbitState GetCartesianGCRFStateAtEpoch(Real epoch) override;

    VecX GetStateVec() {
      Vec6 rv = rv_->GetVec();
      Vec2 clk = GetClockState().GetVec();
      VecX state(8);
      state << rv, clk;
      return state;
    }
  };

}  // namespace lupnt
