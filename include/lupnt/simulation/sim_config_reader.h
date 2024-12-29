/**
 * @file sim_config_reader.h
 * @author Stanford NAV LAB
 * @brief  Simulation configuration file reader
 * @version 0.1
 * @date 2024-12-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include <yaml-cpp/yaml.h>
#include "lupnt/core/definitions.h"
#include "lupnt/agents/agent.h"
#include "lupnt/agents/spacecraft.h"
#include "lupnt/agents/gnss_constellation.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/measurements/space_channel.h"
#include "lupnt/measurements/gnss_channel.h"
#include "lupnt/physics/clock.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/numerics/filters.h"
#include "lupnt/apps/state_estimation_app.h"

namespace lupnt {

    struct TimeConfig {
        double epoch0 = 0.0;
        double dt_integ = 0.0;
        double dt_meas = 0.0;
        double tf = 0.0;
        double n_orbit = 0;
    };

    class ConfigReader {

        protected:
            TimeConfig time_config_;
            Ptr<GnssConstellation> gnss_;
            std::vector<Ptr<Spacecraft>> satellites_;
            std::map<std::string, Ptr<Spacecraft>> satellites_map_;
            std::map<std::string, Ptr<GnssChannel>> gnss_channels_map_;
            std::map<std::string, Ptr<SpaceChannel>> space_channels_map_;
            std::map<std::string, Ptr<IDynamics>> dynamics_map_;
            std::map<std::string, Ptr<ClockDynamics>> clock_dynamics_map_;
            std::map<std::string, Ptr<IState>> state_map_;
            std::map<std::string, Ptr<IFilter>> filter_map_;

            // Yaml node loaders
            void LoadTimeConfig(YAML::Node time_node, bool print_val=false);
            void LoadGnssConfig(YAML::Node gnss_node, bool print_val=false);
            void LoadChannelConfig(YAML::Node channel_node, bool print_val=false);
            void LoadDynamicsConfig(YAML::Node dynamics_node, bool print_val=false);
            void LoadSatelliteConfig(YAML::Node satellite_node, bool print_val=false);
            void LoadFiltersConfig(YAML::Node filter_node, bool print_val=false);

            // state, device, dynamics generators
            void CreateSatCommDevice(YAML::Node device_node, std::string device_name, Ptr<Spacecraft> sat);
            template<typename T> Ptr<NBodyDynamics<T>> CreateNBodyDynamics(const YAML::Node child_node, std::string name, bool print_val);
            void CreateSatState(YAML::Node state_node, std::string state_name, Ptr<Spacecraft> sat, 
                                bool &orbit_defined, bool &clock_defined);
            void CreateOrbitState(YAML::Node state_node, Ptr<Spacecraft> sat, std::string prior_field);
            void CreateClockState(YAML::Node state_node, Ptr<Spacecraft> sat, std::string prior_field);
            Ptr<JointState> CreateJointState(YAML::Node state_node, std::string field_prior);

            // util functions
            template<typename T> T LoadRequiredField(const YAML::Node& node, const std::string& field_name);
            template<>  Real LoadRequiredField(const YAML::Node& node, const std::string& field_name);
            template<typename T> Ptr<T> FindMap(std::map<std::string, Ptr<T>> &map, std::string name);
            Frame FindStrFrameMap(const std::map<std::string, Frame> &map, std::string name);
            NaifId FindFrameCenterMap(const std::map<Frame, NaifId> &map, Frame name);
            ClockModel FindClockModel(std::string key);
            template<typename T> bool InsertMap(std::map<std::string, Ptr<T>> &map, std::string name, Ptr<T> obj);

        public:
            ConfigReader() = default;

            ConfigReader(std::string filename, bool print_val=false) {
                LoadConfigYaml(filename, print_val);
            }

            ~ConfigReader() = default;

            void LoadConfigYaml(std::string filename, bool print_val=false);

            // Getters
            TimeConfig GetTimeConfig() { return time_config_; }
            Ptr<GnssConstellation> GetGnssConstellation() { return gnss_; }
            std::map<std::string, Ptr<Spacecraft>> GetSatellites() { return satellites_map_; }
            std::map<std::string, Ptr<GnssChannel>> GetGnssChannels() { return gnss_channels_map_; }
    };
} // namespace lupnt

