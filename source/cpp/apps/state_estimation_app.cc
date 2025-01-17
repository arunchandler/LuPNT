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

#include "lupnt/apps/state_estimation_app.h"

namespace lupnt {

  // without any parameters
  void JointState::PushBackStateAndDynamics(Ptr<IState> state, Ptr<IDynamics> dynamics,
                                            Ptr<FilterProcessNoiseFunction> proc_noise_func) {
    state_vec_.push_back(state);
    state_vec_size_ += state->GetSize();
    state_types_ += 1;
    state_sizes_.push_back(state->GetSize());

    Ptr<DynamicsWithParams> dyn_params = MakePtr<DynamicsWithParams>(dynamics);
    dynamics_vec_.push_back(dyn_params);

    // update the internal state vector
    state_vec_value_.resize(state_vec_size_);
    int cur_idx = 0;
    for (int i = 0; i < state_types_; i++) {
      for (int j = 0; j < state_vec_[i]->GetSize(); j++) {
        state_vec_value_(cur_idx) = state_vec_[i]->GetValue(j);
        cur_idx++;
      }
    }

    // Add the process noise function
    proc_noise_vec_.push_back(proc_noise_func);
  }

  // Todo: Pushback also parameters and proc noise unctions
  void JointState::PushBackStateAndDynamics(Ptr<IState> state, Ptr<DynamicsWithParams> dynamics,
                                            Ptr<FilterProcessNoiseFunction> proc_noise_func) {
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

    // Add the process noise function
    proc_noise_vec_.push_back(proc_noise_func);
  };

  FilterDynamicsFunction JointState::GetFilterDynamicsFunction() {
    FilterDynamicsFunction dynfunc = [this](VecX x, Real t_curr, Real t_end, MatXd* Phi = nullptr) {
      std::vector<Ptr<IState>> state_vec = this->GetJointState();

      // Iterate for each dynamics and corresponding state (e.g. orbit and
      // dynamics)
      int start_idx = 0;
      int state_vec_n = state_vec.size();
      int dynamics_size = this->dynamics_vec_.size();

      int state_size = 0;

      if (Phi != nullptr) {
        Phi->resize(state_vec_size_, state_vec_size_);
        Phi->setZero();
      }

      for (int i = 0; i < dynamics_size; i++) {
        state_size = this->state_sizes_[i];
        VecX x_seg(state_size);
        VecX x_seg_next(state_size);
        for (int j = 0; j < state_size; j++) {
          x_seg(j) = x(start_idx + j);
        }

        if (Phi != nullptr) {
          MatXd Phi_tmp(state_size, state_size);
          x_seg_next = this->dynamics_vec_[i]->Propagate(x_seg, t_curr, t_end, &Phi_tmp);
          Phi->block(start_idx, start_idx, state_size, state_size) = Phi_tmp;
        } else {
          // Jacobian not needed (e.g. UKF)
          x_seg_next = this->dynamics_vec_[i]->Propagate(x_seg, t_curr, t_end, nullptr);
        }

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

  FilterProcessNoiseFunction JointState::GetFilterProcessNoiseFunction() {
    FilterProcessNoiseFunction proc_noise_func = [this](VecX x, Real t_curr, Real t_end) {
      MatXd Q = MatXd::Zero(state_vec_size_, state_vec_size_);
      int start_idx = 0;
      int state_size = 0;
      int dynamics_size = this->dynamics_vec_.size();

      for (int i = 0; i < dynamics_size; i++) {
        state_size = this->state_sizes_[i];
        MatXd Q_tmp = (*proc_noise_vec_[i])(x.segment(start_idx, state_size), t_curr, t_end);
        Q.block(start_idx, start_idx, state_size, state_size) = Q_tmp;
        start_idx += state_size;
      }

      return Q;
    };

    return proc_noise_func;
  };
};  // namespace lupnt
