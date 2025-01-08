/**
 * @file dynamics_surface.cc
 * @author Stanford NAV Lab
 * @brief  Dynamics of objects on Surface (Rover, GroundStation)
 * @version 0.1
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/physics/body.h"

namespace lupnt {

  SurfaceStaticDynamics::SurfaceStaticDynamics(const NaifId body_id,
                                               const Frame state_frame = Frame::NONE)
      : body_id_(body_id), state_frame_(state_frame) {
    // Set the dynamics frame to the body-fixed frame
    dynamics_frame_ = GetBodyFixedFrameName(body_id);

    if (state_frame_ == Frame::NONE) {
      // Set the state frame to the body-fixed frame
      state_frame_ = dynamics_frame_;
    }
  }

  Vec6 SurfaceStaticDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    // Define lambda function
    Frame state_frame = state_frame_;
    Frame dynamics_frame = dynamics_frame_;

    auto propfunc = [t0, tf, state_frame, dynamics_frame](const Vec6 &x) {
      if (state_frame == Frame::NONE || dynamics_frame == Frame::NONE) {
        throw std::runtime_error("Frame not set for SurfaceStaticDynamics");
      }

      if (state_frame == dynamics_frame) {
        return x;  // No conversion needed
      }

      // Convert the state from GCRF to dynamics frame
      Vec6 x_bf = ConvertFrame(t0, x, state_frame, dynamics_frame);

      // The dynamics is stationary

      // Convert the state back to GCRF frame
      Vec6 xf_out = ConvertFrame(tf, x_bf, dynamics_frame, state_frame);

      return xf_out;
    };

    // Get stm
    Vec6 xf;
    if (stm != nullptr) {
      Vec6 x0_tmp = x0.cast<double>();
      *stm = jacobian(propfunc, wrt(x0_tmp), at(x0_tmp), xf);
    } else {
      xf = propfunc(x0);
    }

    return xf;
  }

  OrbitState SurfaceStaticDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                   Mat6d *stm) {
    Vec6 x0 = state.GetVec();
    Vec6 xf = Propagate(x0, t0, tf, stm);
    OrbitState new_state = state.CreateCopyWithValue(xf);
    return new_state;
  }

}  // namespace lupnt
