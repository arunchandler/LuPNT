/**
 * @file measurement_simulation.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-12-24
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/scheduler.h"
#include "lupnt/numerics/filters.h"

namespace lupnt {

  class Simulation {
  private:
    std::string sim_name_ = "Base";
    std::vector<Ptr<Agent>> agents_;
    Scheduler scheduler_;

  public:
    Simulation() = default;
    ~Simulation() = default;

    void AddAgent(Ptr<Agent> agent) { agents_.push_back(agent); }
    void RunSimulation(Real end_time) { scheduler_.RunSimulation(end_time); }

    // Print
    virtual void UpdateStates(Real epoch) = 0;
    virtual void SetupEnvironment(std::string filename) = 0;
    virtual void PrintStatus() = 0;
  };

  class NavSimulation : public Simulation {
  private:
    std::string sim_name_ = "MeasurementBase";
    std::vector<Ptr<Agent>> agents_;
    std::vector<IFilter> filters_;
    Scheduler scheduler_;

  public:
    NavSimulation() = default;
    ~NavSimulation() = default;

    void UpdateStates(Real epoch) override;
    virtual void SetupEnvironment(std::string filename) override;
    virtual void PrintStatus() override;
  };

};  // namespace lupnt