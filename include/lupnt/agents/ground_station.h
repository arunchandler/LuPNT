/**
 * @file ground_stations.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include <lupnt/agents/agent.h>

namespace lupnt {

  class GroundStation : public Agent {

  protected:
    Real latitude_;
    Real longitude_;

  public:
    GroundStation(NaifId body_id, Real lat, Real lon) : Agent() { 
        SetIsBodyFixed(true); 
        SetBodyId(body_id);
        // Create body for corresponding body_id
        BodyData body_data = GetBodyData(body_id);

        SetDynamics(std::make_shared<SurfaceStaticDynamics>(body_id, body_data.fixed_frame));
        latitude_ = lat;
        longitude_ = lon;
    };

  private:
    Vec3d pos_;
  };
} // namespace lupnt
