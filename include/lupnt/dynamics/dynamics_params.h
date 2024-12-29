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
    enum DynamicsParamsOption {
        TrueFixed,  // Use the true value 
        Estimated,  // Estimate the value
        Consider    // Consider the value uncertainty 
    };

    /**
     * @brief Dynamics with parameters
     * 
     */
    class DynamicsWithParams : public IDynamics {
        private:
            bool use_params_ = false;
            Ptr<IDynamics> dynamics_;
            DynamicsParam params_;

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

            DynamicsWithParams(Ptr<IDynamics> dyn,
                               DynamicsParam params) {
                dynamics_ = dyn;
                use_params_ = true;
                params_ = params;
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
            Ptr<IState> PropagateState(const Ptr<IState>& state, Real t0, Real tf, MatXd* stm = nullptr) override {
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