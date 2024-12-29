/**
 * @file dynamics_params.cc
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-12-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "lupnt/dynamics/dynamics_params.h"

namespace lupnt {

    VecX DynamicsWithParams::Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm) {

        // To compute the jacobians, x_with_params need to capture all estimated and considered parameters
        auto func = [=](const VecX& x_with_est_params) {

            VecX x = x_with_est_params.head(x0.size() - params_.size());
            
            if (use_params_) {
                // extract the parameters from x_with_est_params
                dynamics_->SetParams(params_);
            }
            VecX xf = dynamics_->Propagate(x, t0, tf, nullptr); // Do not compute stm here
            
            if (use_params_) {
                params_ = dynamics_->GetParams(); // Get the updated parameters (often its constant)
                // Todo: Set params to the last elements of x_with_est_params
            }
            return xf;
        };

        VecX x0_tmp = x0.cast<double>();
        VecX xf;
        if (stm != nullptr) {
            *stm = jacobian(func, wrt(x0_tmp), at(x0_tmp), xf);
        } else {
            xf = func(x0_tmp);
        }

        return xf;
    };


}