/**
 * @file dynamics_params.h
 * @author Stanford NAV LAB
 * @brief Interface for dynamics parameters
 * @version 0.1
 * @date 2024-12-28
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/dynamics/dynamics.h"

namespace lupnt {

  // ****************************************************************************
  // Dynamics Parameters Interface
  // ****************************************************************************

  /**
   * @brief Parameter option for the filters
   *
   */
  enum ParamsEstOption {
    TrueFixed,  // Use the true value
    Estimated,  // Estimate the value
    Consider    // Consider the value uncertainty
  };

  extern const std::map<ParamsEstOption, std::string> params_est_option2string;
  extern const std::map<std::string, ParamsEstOption> string2params_est_option;

  /**
   * @brief Dynamics parameters with parameter estimation option
   *
   */
  typedef std::map<std::string, std::pair<ParamsEstOption, VecX>> DynamicsParamWithOptionMap;

  class DynamicsParamWithOption {
  private:
    DynamicsParamWithOptionMap params_;

  public:
    /**
     * @brief Construct a new Dynamics Param With Option object, Set option tp TrueFixed
     *
     * @param dyn
     */
    DynamicsParamWithOption() {}  // Empty constructor

    DynamicsParamWithOption(DynamicsParam params) {
      for (auto const& [key, val] : params) {
        params_[key] = std::make_pair(ParamsEstOption::TrueFixed, val);  // Set default to TrueFixed
      }
      std::cout << "[DynamicsParamWithOption] Parameters are set to TrueFixed by default"
                << std::endl;
    }

    DynamicsParamWithOptionMap GetMap() { return params_; }

    void AddParam(std::string key, ParamsEstOption option, VecX value) {
      params_[key] = std::make_pair(option, value);
    }

    void UpdateParamValue(std::string key, VecX value) {
      if (params_.find(key) != params_.end()) {
        params_[key].second = value;
      }
    }

    void SetOption(std::string key, ParamsEstOption option) {
      if (params_.find(key) != params_.end()) {
        params_[key].first = option;
      }
    }

    void PrintKeys() {
      for (auto const& [key, val] : params_) {
        std::cout << key << ", ";
      }
      std::cout << std::endl;
    }

    DynamicsParam GetDynamicsParam() {
      DynamicsParam params;
      for (auto const& [key, val] : params_) {
        params[key] = val.second;
      }
      return params;
    }
  };

  /**
   * @brief Dynamics with parameters
   *
   */
  class DynamicsWithParams : public IDynamics {
  protected:
    bool use_params_ = false;
    Ptr<IDynamics> dynamics_;
    DynamicsParamWithOption params_;

  public:
    /**
     * @brief Construct a new Dynamics With dynamics only -> no parameters
     *
     * @param dyn
     */
    DynamicsWithParams(Ptr<IDynamics> dyn) {
      use_params_ = false;
      dynamics_ = dyn;
    }

    DynamicsWithParams(Ptr<IDynamics> dyn, DynamicsParam params) {
      dynamics_ = dyn;
      use_params_ = true;
      params_ = DynamicsParamWithOption(params);
    }

    DynamicsWithParams(Ptr<IDynamics> dyn, DynamicsParamWithOption params) {
      dynamics_ = dyn;
      use_params_ = true;
      params_ = params;
    }

    DynamicsWithParams(Ptr<IDynamics> dyn, DynamicsParam params, ParamsEstOption option) {
      dynamics_ = dyn;
      use_params_ = true;
      params_ = DynamicsParamWithOption(params);
      for (auto const& [key, val] : params) {
        params_.SetOption(key, option);  // set all to the same option
      }
    }

    void SetParamEstOption(std::string key, ParamsEstOption option) {
      params_.SetOption(key, option);
    }

    Ptr<IDynamics> GetDynamics() { return dynamics_; }

    DynamicsParam GetParams() { return params_.GetDynamicsParam(); }

    int GetTotalParamSize() {
      int total_size = 0;
      for (auto const& [key, val] : params_.GetMap()) {
        total_size += val.second.size();
      }
      return total_size;
    }

    int GetEstParamSize() {
      int total_size = 0;
      for (auto const& [key, val] : params_.GetMap()) {
        if (val.first == ParamsEstOption::Estimated) {
          total_size += val.second.size();
        }
      }
      return total_size;
    }

    int GetConsiderParamSize() {
      int total_size = 0;
      for (auto const& [key, val] : params_.GetMap()) {
        if (val.first == ParamsEstOption::Consider) {
          total_size += val.second.size();
        }
      }
      return total_size;
    }

    void UpdateParamValues(DynamicsParam params) {
      for (auto const& [key, val] : params) {
        params_.UpdateParamValue(key, val);
      }
    }

    VecX GetEstConsiderParamVector(DynamicsParam dyn_param) {
      VecX param_vec(GetEstParamSize() + GetConsiderParamSize());
      int idx = 0;
      for (auto const& [key, val] : params_.GetMap()) {
        if (val.first == ParamsEstOption::Estimated) {
          param_vec.segment(idx, val.second.size()) = dyn_param[key];
          idx += val.second.size();
        } else if (val.first == ParamsEstOption::Consider) {
          param_vec.segment(idx, val.second.size()) = dyn_param[key];
          idx += val.second.size();
        }
      }
      return param_vec;
    }

    /*
     * @brief Propagate the state
     *
     * @param state
     * @param t0
     * @param tf
     * @param stm
     * @return Ptr<IState>
     */
    Ptr<IState> PropagateState(const Ptr<IState>& state, Real t0, Real tf,
                               MatXd* stm = nullptr) override {
      // warning
      std::cout << "Warning: PropagateState does not consider parameters in stm when using "
                   "DynamicsWithParams"
                << std::endl;
      dynamics_->SetParams(params_.GetDynamicsParam());
      return dynamics_->PropagateState(state, t0, tf, stm);
    }

    /**
     * @brief A function to propagate the dynamics
     *
     * @param x0  concatenated state and parameters (with Either Estimated or Considered)
     * @param t0  initial time
     * @param tf  final time
     * @param stm   state transition matrix
     * @return VecX
     */
    VecX Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm = nullptr) override;
  };

}  // namespace lupnt