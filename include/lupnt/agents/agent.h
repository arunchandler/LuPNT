/**
 * @file agent.h
 * @author Stanford NAV LAB
 * @brief List of agents
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

// C++ includes
#include <memory>
#include <typeinfo>

// lupnt includes
#include "lupnt/core/constants.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/physics/attitude_state.h"
#include "lupnt/physics/clock.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  class ICommDevice;

  /**
   * @brief Agent base class
   *
   */
  class Agent {
  protected:
    static int id_counter_;
    const int id_;
    std::string name_;
    bool is_bodyfixed_;

    NaifId bodyId_;
    Real epoch_;
    Ptr<IState> rv_;
    Ptr<IDynamics> dynamics_;
    Frame dynamics_frame_ = Frame::NONE;

    Ptr<AttitudeState> attitude_;
    std::vector<Ptr<ICommDevice>> devices_;

    ClockState clock_;
    std::unique_ptr<ClockDynamics> clock_dynamics_;

  public:
    Agent() : id_(id_counter_++), clock_(ClockState(2)) {};
    virtual ~Agent() = default;

    // Getters
    Real GetEpoch() const { return epoch_; }
    std::string GetName() const { return name_; }
    NaifId GetBodyId() const { return bodyId_; }
    bool IsBodyFixed() const { return is_bodyfixed_; }
    Ptr<IState> GetRvState() const { return rv_; }
    Ptr<IDynamics> GetDynamics() const { return dynamics_; }
    Ptr<AttitudeState> GetAttitudeState() const { return attitude_; }
    ClockState GetClockState() const { return clock_; }

    // Setters
    void SetIsBodyFixed(bool is_bodyfixed) { is_bodyfixed_ = is_bodyfixed; }
    void SetDynamicsFrame(Frame frame) { dynamics_frame_ = frame; }
    void GetDynamicsFrame(Frame frame) { dynamics_frame_ = frame; }
    void SetRvState(Ptr<IState> rv) { rv_ = rv; }
    void SetDynamics(Ptr<IDynamics> dyn) { dynamics_ = dyn; }
    void SetEpoch(Real epoch) { epoch_ = epoch; }
    void SetBodyId(NaifId bodyId) { bodyId_ = bodyId; }
    void SetClock(ClockState clk) { clock_ = clk; }
    void SetClockDynamics(ClockDynamics& clock_dyn) {
      clock_dynamics_ = std::make_unique<ClockDynamics>(clock_dyn);
    }

    // Comm Device
    void AddDevice(Ptr<ICommDevice> device) { devices_.push_back(device); }

    Ptr<Transmitter> GetTransmitter() {
      for (auto device : devices_) {
        if (device->txrx == "tx") {
          return std::dynamic_pointer_cast<Transmitter>(device);
        }
      }
      return nullptr;
    }

    Ptr<Receiver> GetReceiver() {
      for (auto device : devices_) {
        if (device->txrx == "rx") {
          return std::dynamic_pointer_cast<Receiver>(device);
        }
      }
      return nullptr;
    }

    Ptr<Transponder> GetTransponder() {
      for (auto device : devices_) {
        if (device->txrx == "txrx") {
          return std::dynamic_pointer_cast<Transponder>(device);
        }
      }
      return nullptr;
    }

    // Cartesian OrbitState at epoch in GCRF frame
    virtual CartesianOrbitState GetCartesianGCRFStateAtEpoch(Real epoch) = 0;

    /**
     * @brief Propagate the agent's state to the given epoch
     *
     * @param epoch  Epoch to propagate to
     */
    void Propagate(const Real epoch);

    /**
     * @brief Get the Rv State at epoch, wihtout changing the agent's epoch and
     * state
     *
     * @param epoch  Epoch to get the state at
     * @return VecX
     */
    VecX GetRvStateAtEpoch(const Real epoch);

    /**
     * @brief Propagate the state to the given epoch
     *
     * @param epoch0    Initial epoch
     * @param x         Initial state
     * @param epoch     Final epoch
     * @param inout_frame   Frame of the input state
     * @return
     */
    VecX PropagateRvState(const Real epoch0, const Vec6& x0, const Real epoch, Frame inout_frame);

    /**
     * @brief Get the Clock State at epoch, wihtout changing the agent's epoch and
     * state
     *
     * @param epoch     Epoch to get the state at
     * @param with_noise    Propagate with noise
     * @return ClockState
     */
    ClockState GetClockStateAtEpoch(const Real epoch, bool with_noise = true);

    /**
     * @brief Propagate the state to the given epoch
     *
     * @param epoch0    Initial epoch
     * @param x         Initial state
     * @param epoch     Final epoch
     * @return
     */
    VecX PropagateClockState(const Real epoch0, const Vec2& x0, const Real epoch);

    /**
     * @brief Get the State Vec object
     *
     * @return VecX
     */
    VecX GetClockStateVecAtEpoch(const Real epoch, bool with_noise = true);
  };

};  // namespace lupnt
