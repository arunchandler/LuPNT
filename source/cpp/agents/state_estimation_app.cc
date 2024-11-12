/**
 * @file state_estimation_app.cc
 * @author Stanford NAV Lab
 * @brief  Base class for State Estimation Application
 * @version 0.1
 * @date 2024-11-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "lupnt/agents/state_estimation_app.h"

namespace lupnt {

    void JointState::PushBackStateAndDynamics(IState* state, IDynamics* dynamics) {
      state_vec_.push_back(state);
      state_vec_size_ += state->GetSize();
      state_types_ += 1;
      state_sizes_.push_back(state->GetSize());

      dynamics_vec_.push_back(dynamics);

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

    FilterDynamicsFunction JointState::GetFilterDynamicsFunction() {
      FilterDynamicsFunction dynfunc = [this](VecX x, Real t_curr, Real t_end, MatXd& Phi) {
        std::vector<IState*> state_vec = this->GetJointState();
        Phi.resize(state_vec_size_, state_vec_size_);
        Phi.setZero();

        // Iterate for each dynamics and corresponding state (e.g. orbit and
        // dynamics)
        int start_idx = 0;
        int state_vec_n = state_vec.size();
        int dynamics_size = this->dynamics_vec_.size();

        int state_size = 0;

        for (int i = 0; i < dynamics_size; i++) {
          state_size = this->state_sizes_[i];
          MatXd Phi_tmp(state_size, state_size);
          VecX x_seg(state_size);
          VecX x_seg_next(state_size);
          for (int j = 0; j < state_size; j++) {
            x_seg(j) = x(start_idx + j);
          }
          x_seg_next = this->dynamics_vec_[i]->Propagate(x_seg, t_curr, t_end, &Phi_tmp);
          Phi.block(start_idx, start_idx, state_size, state_size) = Phi_tmp;
          for (int j = 0; j < state_size; j++) {
            x(start_idx + j) = x_seg_next(j);
          }
          // Add states
          start_idx += state_size;
        }

        return x;
      };

      return dynfunc;
    };
}; // namespace lupnt
