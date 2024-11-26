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
} // namespace lupnt
