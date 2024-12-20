/**
 * @file spacecraft.cc
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/agents/spacecraft.h"

namespace lupnt {
  /********************** SpaceCraft  ********************************/

  CartesianOrbitState Spacecraft::GetCartesianGCRFStateAtEpoch(Real epoch) {
    Ptr<OrbitState> state = GetOrbitState();

    Real current_epoch = GetEpoch();
    Ptr<NumericalOrbitDynamics> dynamics
        = std::dynamic_pointer_cast<NumericalOrbitDynamics>(GetDynamics());

    Real GM = GetBodyData(GetBodyId()).GM;

    if (epoch != current_epoch) {
      // set dt
      OrbitState prop_state = dynamics->PropagateState(*state, current_epoch, epoch);
      // Create a pointer to new_state
      Ptr<OrbitState> new_state = MakePtr<OrbitState>(prop_state);

      Ptr<CartesianOrbitState> cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
          ConvertOrbitStateRepresentation(new_state, OrbitStateRepres::CARTESIAN, GM));
      return ConvertOrbitStateFrame(*cartOrbitState, epoch, Frame::GCRF);
    } else {
      // No need to propagate
      Ptr<CartesianOrbitState> cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
          ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, GM));
      return ConvertOrbitStateFrame(*cartOrbitState, epoch, Frame::GCRF);
    }
  }
}  // namespace lupnt
