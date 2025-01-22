/**
 * @file keplarian_ephemeris.h
 * @author Stanford NAV LAB
 * @brief Ephemeris based on keplarian state
 * @version 0.1
 * @date 2025-01-12
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/physics/frame_converter.h"
#include <map>

namespace lupnt {

    struct EphemFitOption {
        int max_iter = 10;
        double lm_lambda0 = 1e-3;
        double lm_down = 9.0;
        double lm_up = 11.0;    
        double eps_cost = 1e-11;
        double eps_lm = 1e-1;
        double lm_max = 1e7;
        double lm_min = 1e-7;
        bool debug = false;
    };

    struct EphemFitResult {
        VecXd t;
        VecXd x;
        VecXd x0;
        MatXd pos_true;
        MatXd vel_true;
        MatXd vel_true_rot;  // Earth rotation included
        MatXd pos_est;
        MatXd vel_est;
        MatXd pos_init;
        MatXd vel_init;
        double cost;
        bool success;
    };

    /**
     * @brief Interface for keplarian ephemeris
     * 
     */
    class IEphemeris {
        protected:
            // State vector
            int ephem_size_ = 0;  // Size of the ephemeris
            std::string name_; // Name of the ephemeris
            std::filesystem::path save_path_; // Save path
            std::vector<std::string> ephemeris_names_;  // State vector names
            NaifId body_id_ = NaifId::MOON;  // Body ID
            bool save_path_set_ = false;

        public:
            IEphemeris() = default;

            virtual ~IEphemeris() = default;

            /**
             * @brief Convert the ephemris to the fixed frame position
             * 
             * @param epoch  Epoch in tai
             * @param ephemeris  Ephemeris vector
             * @return Vec3 Position in the fixed frame
             */
            virtual Vec3 EphemerisToFixedFramePos(Real epoch, VecX ephemeris) = 0;

            /**
             * @brief Convert the ephemeris to fixed frame position and velocity
             * 
             * @param epoch  Epoch in tai
             * @param ephemeris  Ephemeris vector
             * @return Vec6 Position and velocity in the fixed frame
             */
            virtual Vec6 EphemerisToFixedFramePosVel(Real epoch, VecX ephemeris) = 0;

            /**
             * @brief Convert the epheemris to the fixed frame state
             * 
             * @param epoch 
             * @param ephemeris 
             * @param H 
             * @return VecX 
             */
            Vec3 EphemerisToFixedFramePos(Real epoch, VecX ephemeris, MatXd* H);

            /**
             * @brief Get the Initial of the ephemeris 
             * 
             * @param t_fit:  time vector in tai
             * @param fit_arc_bf:  MatX (lent, 6) in body fixed frame
             * @return VecX   Initial guess of the ephemeris
             */
            virtual VecX GetInitialGuess(VecX t_fit, MatX fit_arc_bf) = 0;

            /**
             * @brief  Fit the ephemeris from the ephemeris arc using non-linear least squares
             * 
             * @param t_fit  Time vector (Epoch in tai)
             * @param fit_arc  Arc of the ephemeris (n x 6)
             * @param fit_option  Fit options
             */
            EphemFitResult FitEphemeris(VecX t_fit, MatX fit_arc, Frame frame, EphemFitOption fit_option, bool optimize = true);

            /**
             * @brief Plot the fitted orbit
             * 
             * @param fit_result  Fit result
             * @param is_plot  Plot flag
             * @param filename  Filename for the figure
             * @param plot_init  Plot the initial orbit
             */
            Vec6d EvalFitError(EphemFitResult fit_result, bool is_plot, std::string fig_name, bool plot_init);

            /**
             * @brief Get the ephemeris map from the ephemeris vector
             * 
             * @param ephemeris 
             * @return std::map<std::string, Real> 
             */
            std::map<std::string, Real> GetEphemerisMap(VecX ephemeris) const { 
                std::map<std::string, Real> ephemeris_map;
                for (int i = 0; i < ephem_size_; i++) {
                    ephemeris_map[ephemeris_names_[i]] = ephemeris(i);
                }
                return ephemeris_map; 
            };

            int GetEphemIndex(std::string field) {
                // find the index for field
                for (int i = 0; i < ephem_size_; i++) {
                    if (ephemeris_names_[i] == field) {
                        return i;
                    }
                }
                std::cerr << "Field: " << field << " not found in ephemeris" << std::endl;
                return 1000;
            };

            std::vector<std::string> GetEphemerisNames() const { return ephemeris_names_; }
            int GetEphemerisSize() const { return ephem_size_; }
            NaifId GetBodyId() const { return body_id_; }

            std::filesystem::path GetSavePath() const { return save_path_; }
            void SetSavePath(std::filesystem::path path) { 
                save_path_ = path;
                save_path_set_ = true; 
            }
    };

    /**
     * @brief Ephemeris for GPS satellites
     * 
     */
    enum class KepEphemType {
        GPS16,
        EPH18,
        EPH20,
        EPH22,
    };

    class KeplarianEphemeris : public IEphemeris {

        private:
            bool use_dotu_ = false;
            bool use_dotr_ = false;
            bool use_cl_ = false;
            bool use_cr2_ = false;

        public:
            KeplarianEphemeris(NaifId body_id, KepEphemType type);

            Vec3 EphemerisToFixedFramePos(Real epoch, VecX ephemeris) override;
            Vec6 EphemerisToFixedFramePosVel(Real epoch, VecX ephemeris) override;

            VecX GetInitialGuess(VecX t_fit, MatX fit_arc_bf) override;

            /**
             * @brief Compute r_k, u_k, i_k, lambda_k from the ephemeris
             * 
             * @param t  Time
             * @param ephemeris  Ephemeris vector
             * @return Vec4   r_k, u_k, i_k, lambda_k
             */
            std::map<std::string, Real> ComputeInternalParams(Real t, VecX ephemeris, bool compute_vel = false);

            /**
             * @brief Get the Initial Guess for Cr and Cs coefficients
             * 
             * @param fit_t   Time vector
             * @param fit_arc  Arc of the ephemeris
             * @param t_ref   Reference time
             * @param coe_ref  keplarian elements at ref time
             * @param frame   Frame of the ephemeris
             * @return VecX   Cuc, Cus, Crc, Crs, delta_n
             */
            std::map<std::string, Real> ComputeInitHarmCoeff(VecX fit_t, MatX fit_arc_bf, Real t_ref, Vec6 coe_ref, bool plot_result = false);


    };
}


