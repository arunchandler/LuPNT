/**
 * @file StateEstimationApp.h
 * @author Stanford NAV LAB
 * @brief Base class for State Estimation Application
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <memory>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/application.h"
#include "lupnt/core/file.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/gnss_measurement.h"
#include "lupnt/measurements/gnss_receiver.h"
#include "lupnt/numerics/filters.h"
#include "lupnt/physics/state.h"

namespace lupnt {

  /**
   * @brief Stack of multiple state types (example: orbit and clock)
   *
   */
  class JointState {
  private:
    std::vector<IState*> state_vec_;
    std::vector<IDynamics*> dynamics_vec_;
    std::vector<std::vector<int>> dynamics_to_state_map_;

    VecX state_vec_value_;
    std::vector<int> state_sizes_;
    int state_vec_size_ = 0;
    int state_types_ = 0;

  public:
    JointState(){
      state_vec_value_ = VecX::Zero(0);
      state_vec_size_ = 0;
      state_types_ = 0;
      state_sizes_.clear();
    };

    JointState(std::vector<IState*> state_vec) {
      int state_vec_size = 0;
      for (int i = 0; state_vec.size(); i++) {
        state_vec_.push_back(state_vec[i]);
        state_vec_size += state_vec_[i]->GetSize();
        state_types_ += 1;
        state_sizes_.push_back(state_vec_[i]->GetSize());
      }
      state_vec_size_ = state_vec_size;

      // update the internal state vector
      state_vec_value_.resize(state_vec_size_);
      int cur_idx = 0;
      for (int i = 0; i < state_types_; i++) {
        for (int j = 0; j < state_vec_[i]->GetSize(); j++) {
          state_vec_value_(cur_idx) = state_vec_[i]->GetValue(j);
          cur_idx++;
        }
      }
    };

    int GetSize() const { return state_vec_size_; };
    std::vector<IState*> GetJointState() { return state_vec_; };
    VecX GetJointStateValue() { return state_vec_value_; };

    void PushBackStateAndDynamics(IState* state, IDynamics* dynamics);

    FilterDynamicsFunction GetFilterDynamicsFunction();

  };

  class StateEstimationApp : public Application {
  private:
    double epoch0_;  // start epoch in TAU
    double epoch_;   // curent epoch
    double t_;       // Current time [s]

    Ptr<Agent> agent_;
    Ptr<IFilter> filter_;
    FilterDynamicsFunction dynamics_func_;
    FilterMeasurementFunction meas_func_;
    JointState state_vec_;

  public:
    void SetAgent(Ptr<Agent> agent) { this->agent_ = agent; }

    double GetInitialEpoch() { return epoch0_; };
    double GetCurrentEpoch() { return epoch_; };
    double GetCurrrentTime() { return t_; };

    void Setup();
    void Step(double t_end);  // execute to step t
  };
};  // namespace lupnt
