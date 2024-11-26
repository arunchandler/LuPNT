/**
 * @file rover.h
 * @author Stanford NAV LAB
 * @brief  Rover Agent
 * @version 0.1
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include <lupnt/agents/agent.h>

namespace lupnt {
  /**
   * @brief Rover Agent
   *
   */
  class Rover : public Agent {
  public:
    Rover() : Agent() { SetIsBodyFixed(true); };
  };

  class GroundStation : public Agent {
  public:
    GroundStation() : Agent() { SetIsBodyFixed(true); };

    void SetPosition(Vec3d pos) { pos_ = pos; }
    Vec3d GetPosition() { return pos_; }

  private:
    Vec3d pos_;
  };
} // namespace lupnt
