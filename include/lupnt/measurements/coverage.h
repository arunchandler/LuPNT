/**
 * @file coverage.h
 * @author Stanford NAV Lab
 * @brief  Coverage calculations
 * @version 0.1
 * @date 2024-12-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {
    
    class SurfaceCoverage {

    private:
        NaifId body_id_;
        double elev_mask_deg_ = 0.0;
        bool computed_target_vec_ = false;

    public:
        VecXd epochs_tai_;
        VecXd lat_vec_;
        VecXd lon_vec_;
        MatXd lat_mesh_;
        MatXd lon_mesh_;
        std::vector<MatXd> az_target_;
        std::vector<MatXd> el_target_;

        SurfaceCoverage(NaifId body_id) : body_id_(body_id) {};

        ~SurfaceCoverage() = default;

        void CreateLatLonMesh(VecXd lat_vec, VecXd lon_vec);

        void SetTarget(VecXd epochs_tai, const VecXd& target_rv);

        void SetTarget(VecXd epochs_tai, const MatXd& target_rv);
        

        /**
         * @brief  Compute the surface coverage of a satellite
         * 
         * @param epochs_tai  List of TAI epochs [s] (N,)
         * @param sat_rv      Satellite position and velocity [km, km/s] (N, 6)
         * @param lat_mesh    Latitude meshgrid [deg] (nlon x nlat)
         * @param lon_mesh    Longitude meshgrid [deg] (nlon x nlat)
         * @param sat_frame   Satellite frame 
         * @param body_id     NAIF ID of the body
         * @param elev_mask_deg   Elevation mask [deg]
         * @return MatXd      Surface coverage [0, 1] (nlon x nlat)
         */
        MatXd ComputeSurfaceCoverage(const double elev_mask_deg);

    };

}  // namespace lupnt
