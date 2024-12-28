/**
 * @file sim_config_reader.cc
 * @brief Implementation of ConfigReader for reading YAML config.
 * @date 2024-12-25
 */

#include "lupnt/simulation/sim_config_reader.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/physics/body.h"
#include "lupnt/agents/agent.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/measurements/gnss_receiver.h"
#include "lupnt/physics/time_converter.h"

namespace lupnt {

    //==============================================================================
    // InsertMap and FindMap
    //==============================================================================

    template<typename T>
    bool ConfigReader::InsertMap(std::map<std::string, Ptr<T>> &container,
                                 std::string key, 
                                 Ptr<T> value) {
        auto [iter, inserted] = container.insert({key, value});
        if (!inserted) {
            std::cerr << "[InsertMap] Key '" << key 
                      << "' already exists.\n";
            throw std::runtime_error("Key '" + key + "' already exists.");
            // return false; // unreachable after throw
        }
        return true;
    }

    template<typename T> 
    Ptr<T> ConfigReader::FindMap(std::map<std::string, Ptr<T>> &container,
                                 std::string key) {
        auto iter = container.find(key);
        if (iter == container.end()) {
            std::cerr << "[FindMap] Key '" << key 
                      << "' not found.\n";
            throw std::runtime_error("Key '" + key + "' not found.");
        }
        return iter->second;
    }

    Frame ConfigReader::FindStrFrameMap(const std::map<std::string, Frame> &container,
                                        std::string key) {
        auto iter = container.find(key);
        if (iter == container.end()) {
            std::cerr << "[FindMap] Key '" << key 
                      << "' not found.\n";
            throw std::runtime_error("Key '" + key + "' not found.");
        }
        return iter->second;
    }

    NaifId ConfigReader::FindFrameCenterMap(const std::map<Frame, NaifId> &container,
                                            Frame key) {
        auto iter = container.find(key);
        if (iter == container.end()) {
            std::cerr << "[FindMap] Key '" << key 
                      << "' not found.\n";
            throw std::runtime_error("Key frame not found.");
        }
        return iter->second;
    }

    ClockModel ConfigReader::FindClockModel(std::string key) {
        if (key == "MicrosemiCsac") {
            return ClockModel::kMicrosemiCsac;
        } else if (key == "Rafs") {
            return ClockModel::kRafs;
        } else if (key == "Uso") {
            return ClockModel::kUso;
        } else if (key == "MiniRafs") {
            return ClockModel::kMiniRafs;
        } else {
            std::cerr << "[FindClockModel] Invalid clock model: " << key << "\n";
            return ClockModel::kUnknown;
        }
    }

    //==============================================================================
    // LoadRequiredField
    //==============================================================================

    template<typename T>
    T ConfigReader::LoadRequiredField(const YAML::Node& node,
                                      const std::string& field_name /*= ""*/)
    {
        if (!node) {
            if (!field_name.empty()) {
                std::cerr << "[LoadRequiredField] Required field \""
                          << field_name << "\" not found.\n";
                throw std::runtime_error("Required field \"" + field_name + "\" not found.");
            } else {
                std::cerr << "[LoadRequiredField] Required node not found.\n";
                throw std::runtime_error("Required node not found.");
            }
        }
        return node.as<T>();
    }

    template<>
    Real ConfigReader::LoadRequiredField(const YAML::Node& node,
                                         const std::string& field_name /*= ""*/)
    {
        double val = LoadRequiredField<double>(node, field_name);
        return Real(val);
    }

    //==============================================================================
    // LoadConfigYaml
    //==============================================================================

    void ConfigReader::LoadConfigYaml(std::string filename, bool print_val) {
        YAML::Node config = YAML::LoadFile(filename);

        // This order is important! 
        LoadTimeConfig(config["time"], print_val);
        LoadDynamicsConfig(config["dynamics"], print_val);
        LoadChannelConfig(config["channels"], print_val);
        LoadGnssConfig(config["gnss"], print_val);
        LoadSatelliteConfig(config["agents"], print_val);
    }

    //==============================================================================
    // LoadTimeConfig
    //==============================================================================

    void ConfigReader::LoadTimeConfig(YAML::Node time_node, bool print_val) {
        // time_node["epoch0"] is an array of 6 elements [year, month, day, hour, min, sec]
        // Convert to TAI epoch in seconds
        std::vector<double> epoch0_vec(6, 0.0);
        YAML::Node epoch0_array = time_node["epoch0"]; 
        for (int i = 0; i < 6; i++) {
            // build an error string like "time.epoch0[0]"
            std::string fieldStr = "time.epoch0[" + std::to_string(i) + "]";
            epoch0_vec[i] = LoadRequiredField<double>(epoch0_array[i], fieldStr);
        }
        Real epoch0_utc = Gregorian2Time(epoch0_vec[0],epoch0_vec[1], epoch0_vec[2],epoch0_vec[3],
                                         epoch0_vec[4], Real(epoch0_vec[5]));
        double epoch0_tai = UTC2TAI(epoch0_utc).val();
        time_config_.epoch0 = epoch0_tai;

        time_config_.dt_integ = LoadRequiredField<double>(time_node["dt_integ"], "time.dt_integ");
        time_config_.dt_meas  = LoadRequiredField<double>(time_node["dt_meas"],  "time.dt_meas");
        time_config_.tf       = LoadRequiredField<double>(time_node["tf"],       "time.tf");
        time_config_.n_orbit  = LoadRequiredField<double>(time_node["n_orbit"],  "time.n_orbit");

        if (print_val) {
            std::cout << "Time Config:\n";
            std::cout << "  Epoch0:    " << time_config_.epoch0 << "\n";
            std::cout << "  dt_integ:  " << time_config_.dt_integ << "\n";
            std::cout << "  dt_meas:   " << time_config_.dt_meas  << "\n";
            std::cout << "  tf:        " << time_config_.tf       << "\n";
            std::cout << "  n_orbit:   " << time_config_.n_orbit  << "\n";
        }

        std::cout << "Set Time Config\n";
    }

    //==============================================================================
    // CreateNBodyDynamics (Templated Helper)
    //==============================================================================

    template<typename T> 
    Ptr<NBodyDynamics<T>> ConfigReader::CreateNBodyDynamics(const YAML::Node child_node, bool print_val) {
        // Use LoadRequiredField for 'integrator'
        std::string integrator = LoadRequiredField<std::string>(
            child_node["integrator"], "dynamics.{}.integrator"
        );

        Ptr<NBodyDynamics<T>> dyn;
        if (integrator == "RKF45") {
            dyn = MakePtr<NBodyDynamics<T>>(IntegratorType::RKF45);
        }
        else if (integrator == "RK4") {
            dyn = MakePtr<NBodyDynamics<T>>(IntegratorType::RK4);
        }
        else {
            std::cerr << "Invalid integrator type: only RK4 or RKF45 accepted.\n";
            throw std::runtime_error("Invalid integrator type for NBodyDynamics.");
        }

        // Bodies
        if (child_node["bodies"]) {
            auto bodiesNode = child_node["bodies"];
            for (std::size_t i = 0; i < bodiesNode.size(); i++) {
                std::string body_name = bodiesNode[i].as<std::string>();

                int n = 0; 
                int m = 0;

                bool sphOrderMatch   = (child_node["sph_order"]
                                        && (child_node["sph_order"].size() == bodiesNode.size()));
                bool sphDegreeMatch  = (child_node["sph_degree"]
                                        && (child_node["sph_degree"].size() == bodiesNode.size()));

                if (sphOrderMatch) {
                    n = child_node["sph_order"][i].as<int>();
                }
                if (sphDegreeMatch) {
                    m = child_node["sph_degree"][i].as<int>();
                }
                // Create BodyT<double> 
                BodyT<T> body = CreateBody<T>(body_name, n, m);
                dyn->AddBody(body);
            }
        }

        // SRP
        if (child_node["srp"] && child_node["srp"]["enable"]) {
            dyn->SetUseSrp(true);

            dyn->SetArea(LoadRequiredField<Real>(child_node["srp"]["area"], "dynamics.{}.srp.area"));
            dyn->SetSrpCoeff(LoadRequiredField<Real>(child_node["srp"]["cr"], "dynamics.{}.srp.cr"));
            dyn->SetMass(LoadRequiredField<Real>(child_node["mass"], "dynamics.{}.mass"));
        }
        // Drag
        if (child_node["drag"] && child_node["drag"]["enable"]) {
            dyn->SetUseDrag(true);

            dyn->SetArea(LoadRequiredField<Real>(child_node["drag"]["area"], "dynamics.{}.drag.area"));
            dyn->SetDragCoeff(LoadRequiredField<Real>(child_node["drag"]["cd"], "dynamics.{}.drag.cd"));
            dyn->SetMass(LoadRequiredField<Real>(child_node["mass"], "dynamics.{}.mass"));
        }

        // Integrator options
        if (child_node["integ_options"]) {
            IntegratorParams iparams;
            iparams.abstol = LoadRequiredField<double>(
                child_node["integ_options"]["abstol"], 
                "dynamics.{}.integ_options.abstol"
            );
            iparams.reltol = LoadRequiredField<double>(
                child_node["integ_options"]["reltol"],
                "dynamics.{}.integ_options.reltol"
            );
            dyn->SetIntegratorParams(iparams);
        }

        return dyn;
    }

    //==============================================================================
    // LoadDynamicsConfig
    //==============================================================================

    void ConfigReader::LoadDynamicsConfig(YAML::Node dynamicsNode, bool print_val) {
        for (YAML::const_iterator it = dynamicsNode.begin(); it != dynamicsNode.end(); ++it) {

            std::string name = it->first.as<std::string>();
            YAML::Node child_node = it->second;

            // We assume 'type' is mandatory
            std::string dyn_type = LoadRequiredField<std::string>(child_node["type"], name + ".type");

            Ptr<IDynamics> dyn;

            if (dyn_type == "NBodyDynamics") {
                // If 'use_real' is mandatory:
                bool use_real = LoadRequiredField<bool>(child_node["use_real"], name + ".use_real");
                if (use_real) {
                    dyn = CreateNBodyDynamics<Real>(child_node, print_val);
                }
                else {
                    dyn = CreateNBodyDynamics<double>(child_node, print_val);
                }
            }
            else if (dyn_type == "ClockDynamics") {
                std::string model_name = LoadRequiredField<std::string>(child_node["model"], 
                                                                        name + ".model");
                dyn = MakePtr<ClockDynamics>(FindClockModel(model_name));
                Ptr<ClockDynamics> clock_dyn_ptr = std::dynamic_pointer_cast<ClockDynamics>(dyn);
                InsertMap<ClockDynamics>(clock_dynamics_map_, name, clock_dyn_ptr);
            }
            else {
                std::cerr << "Invalid Dynamics type at '" << name << "'.\n"
                          << "Currently accepted are NBodyDynamics or ClockDynamics.\n";
            }

            // Insert into map
            InsertMap<IDynamics>(dynamics_map_, name, dyn);
        }

        if (print_val) {
            std::cout << "Dynamics Config:\n";
            for (auto &entry : dynamics_map_) {
                std::cout << "  - Dynamics Name: " << entry.first << "\n";
            }
        }
        std::cout << "Set Dynamics Config\n";
    }

    //==============================================================================
    // LoadChannelConfig
    //==============================================================================

    void ConfigReader::LoadChannelConfig(YAML::Node channelNode, bool print_val) {
        for (YAML::const_iterator it = channelNode.begin(); it != channelNode.end(); ++it) {
            std::string name = it->first.as<std::string>();
            YAML::Node child_node = it->second;

            std::string channel_type = LoadRequiredField<std::string>(
                child_node["type"], 
                name + ".type"
            );

            if (channel_type == "GnssChannel") {
                InsertMap<GnssChannel>(gnss_channels_map_, name, MakePtr<GnssChannel>());
            } else {
                // fallback or other channel
                InsertMap<SpaceChannel>(space_channels_map_, name, MakePtr<SpaceChannel>());
            }
        }
        // If you need print_val for debugging, you can add prints here
    }

    //==============================================================================
    // LoadGnssConfig
    //==============================================================================

    void ConfigReader::LoadGnssConfig(YAML::Node gnssNode, bool print_val) {
        bool is_first_const = true;

        for (YAML::const_iterator it = gnssNode.begin(); it != gnssNode.end(); ++it) {
            std::string name = it->first.as<std::string>();
            YAML::Node child_node = it->second;

            std::string field_prior = "gnss." + name;

            // initialization is mandatory
            std::string init_type = LoadRequiredField<std::string>(
                child_node["initialization"], field_prior + ".initialization"
            );

            if (init_type == "TLE") {
                if (is_first_const) {
                    // Check if we have the required dynamics
                    std::string dyn_key = LoadRequiredField<std::string>(
                        child_node["dynamics"], "gnss.dynamics"
                    );
                    if (dynamics_map_.find(dyn_key) == dynamics_map_.end()) {
                        std::cerr << "Dynamics '" << dyn_key 
                                  << "' not found for GNSS constellation.\n";
                        return;
                    }
                    std::string gnss_type = LoadRequiredField<std::string>(child_node["gnss_type"], field_prior + ".gnss_type");
                    std::string filename  = LoadRequiredField<std::string>(child_node["filename"], field_prior + ".filename");
                    std::string channel   = LoadRequiredField<std::string>(child_node["channel"],  field_prior + ".channel");
                    auto dyn_gnss = FindMap<IDynamics>(dynamics_map_, dyn_key);
                    auto channel_gnss = FindMap<GnssChannel>(gnss_channels_map_, channel);

                    gnss_->InitializeWithTle(gnss_type, filename, dyn_gnss,
                                             channel_gnss, time_config_.epoch0);
                    is_first_const = false;
                } else {
                    // Additional satellites in TLE
                    std::string gnss_type = LoadRequiredField<std::string>(
                        child_node["gnss_type"], field_prior + ".gnss_type"
                    );
                    std::string filename  = LoadRequiredField<std::string>(
                        child_node["filename"], field_prior + ".filename"
                    );
                    gnss_->AddSatellitesWithTle(gnss_type, filename);
                }
            }
            else {
                std::cerr << "Invalid GNSS initialization type: only TLE supported.\n";
            }
        }
    }

    //==============================================================================
    // CreateOrbitState
    //==============================================================================

    void ConfigReader::CreateOrbitState(YAML::Node state_node, Ptr<Spacecraft> sat) {
        // We want to find final frame
        std::string sat_name = sat->GetName();
        std::string prior_key = "satellites." + sat_name + ".states.";

        //  - .frame and .init.frame are required
        std::string final_frame_key = LoadRequiredField<std::string>(
            state_node["frame"], prior_key + "frame"
        );
        Frame frame_final = FindStrFrameMap(string2frame, final_frame_key);

        // .init.frame
        std::string init_frame_key = LoadRequiredField<std::string>(
            state_node["init"]["frame"], prior_key + "init.frame"
        );
        Frame frame_init = FindStrFrameMap(string2frame, init_frame_key);

        // center body from frame_init
        NaifId center_init = FindFrameCenterMap(frame_centers, frame_init);
        NaifId center_final = FindFrameCenterMap(frame_centers, frame_final);

        Real GM = Real(GetBodyGM(center_init));

        // Center Body
        sat->SetBodyId(center_final);

        // Check 'class'
        if (!state_node["class"]) {
            std::cerr << "Class not found for Orbit State.\n";
            return;
        }

        if (state_node["class"].as<std::string>() == "CartesianOrbitState") {
            // Which type is used for 'init'? classical or cartesian?
            if (!state_node["init"]["type"]) {
                std::cerr << "Initialization type not found for Orbit State.\n";
                return;
            }
            std::string init_type = LoadRequiredField<std::string>(state_node["init"]["type"], "orbit.init.type");

            if (init_type == "ClassicalOE") {
                // Build COE
                Vec6 coe_vec = {
                    LoadRequiredField<Real>(state_node["init"]["val"]["a"],        prior_key + "init.val.a"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["e"],        prior_key + "init.val.e"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["i"],        prior_key + "init.val.i") * RAD,
                    LoadRequiredField<Real>(state_node["init"]["val"]["Omega"],    prior_key + "init.val.Omega") * RAD,
                    LoadRequiredField<Real>(state_node["init"]["val"]["w"],        prior_key + "init.val.w") * RAD,
                    LoadRequiredField<Real>(state_node["init"]["val"]["M"],        prior_key + "init.val.M") * RAD
                };
                ClassicalOE coe(coe_vec, frame_init);
                CartesianOrbitState state_init = Classical2Cart(coe, GM);
                CartesianOrbitState final_state = ConvertOrbitStateFrame(state_init, time_config_.epoch0, frame_final);
                sat->SetOrbitState(MakePtr<CartesianOrbitState>(final_state));
            }
            else if (init_type == "CartesianOrbitState") {
                Vec6 state_init_vec = {
                    LoadRequiredField<Real>(state_node["init"]["val"]["rx"], prior_key + "init.val.rx"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["ry"], prior_key + "init.val.ry"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["rz"], prior_key + "init.val.rz"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["vx"], prior_key + "init.val.vx"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["vy"], prior_key + "init.val.vy"),
                    LoadRequiredField<Real>(state_node["init"]["val"]["vz"], prior_key + "init.val.vz")
                };
                CartesianOrbitState state_init(state_init_vec, frame_init);
                CartesianOrbitState final_state =
                    ConvertOrbitStateFrame(state_init, time_config_.epoch0, frame_final);
                sat->SetOrbitState(MakePtr<CartesianOrbitState>(final_state));
            }
            else {
                std::cerr << "Invalid initialization type for Orbit State: must be ClassicalOE or CartesianOrbitState.\n";
            }
        } 
        else {
            std::cerr << "Invalid State class for Orbit State: must be CartesianOrbitState.\n";
        }

        // If there's a 'dynamics' field
        if (state_node["dynamics"]) {
            // we treat as optional or mandatory?
            std::string dyn_key = state_node["dynamics"].as<std::string>();
            auto dyn = FindMap<IDynamics>(dynamics_map_, dyn_key);
            sat->SetDynamics(dyn);
        } else {
            std::cerr << "Dynamics not found for Orbit State.\n";
        }
    }

    //==============================================================================
    // CreateClockState
    //==============================================================================

    void ConfigReader::CreateClockState(YAML::Node state_node, Ptr<Spacecraft> sat) {
        // We expect "val.bias" and "val.drift" at minimum
        Real bias  = Real(LoadRequiredField<double>(state_node["val"]["bias"],  "clock.val.bias"));
        Real drift = Real(LoadRequiredField<double>(state_node["val"]["drift"], "clock.val.drift"));

        VecX clock_vec;
        if (state_node["val"]["drift_rate"]) {
            Real drift_rate = Real(state_node["val"]["drift_rate"].as<double>());
            clock_vec << bias, drift, drift_rate;
        }
        else {
            clock_vec << bias, drift;
        }
        ClockState clock_state(clock_vec);
        sat->SetClock(clock_state);

        // If there's clock dynamics
        if (state_node["dynamics"]) {
            auto clock_dyn = FindMap<ClockDynamics>(clock_dynamics_map_,  state_node["dynamics"].as<std::string>());
            sat->SetClockDynamics(clock_dyn);
        }
    }

    //==============================================================================
    // CreateSatCommDevice
    //==============================================================================

    void ConfigReader::CreateSatCommDevice(YAML::Node device_node, std::string device_name, Ptr<Spacecraft> sat) {
        
        std::string prior_field = "device." + device_name;
        std::string device_class = LoadRequiredField<std::string>(device_node["class"], prior_field + ".class");
        std::string device_channel = LoadRequiredField<std::string>(device_node["channel"], prior_field + ".channel");

        if (device_class == "GnssReceiver") {
            Ptr<GnssReceiver> device;
            std::string rcvr_name = LoadRequiredField<std::string>(device_node["name"], prior_field + ".name");
            std::string attitude_mode = LoadRequiredField<std::string>(device_node["attitude_mode"], 
                                                                        prior_field + ".attitude_mode");
            device = MakePtr<GnssReceiver>(rcvr_name);
            device->SetReceiverAttitudeMode(attitude_mode); 
            // Interface Settings with Agent and Channel                                                
            device->SetChannel(FindMap<GnssChannel>(gnss_channels_map_, device_channel));
            gnss_channels_map_[device_channel]->AddReceiver(device);
            sat->AddDevice(device);
            device->SetAgent(sat);
        }
        else {
            std::cerr << "Invalid Device type for '" << device_name << "': only GnssReceiver supported.\n";
            return;
        }
    }

    //==============================================================================
    // CreateSatState
    //==============================================================================

    void ConfigReader::CreateSatState(YAML::Node state_node, std::string state_name, Ptr<Spacecraft> sat, 
                                      bool &orbit_defined, bool& clock_defined) {
        // "class" is mandatory
        std::string prior_field = "satellites." + sat->GetName() + ".states." + state_name;

        std::string class_name = LoadRequiredField<std::string>(state_node["class"], prior_field + ".class");

        if (class_name == "CartesianOrbitState") {
            if (orbit_defined) {
                std::cerr << "Orbit State already defined for Agent '" 
                            << sat->GetName() << "'\n";
                return;
            }
            CreateOrbitState(state_node, sat);
            orbit_defined = true;

            InsertMap<IState>(state_map_, state_name, sat->GetOrbitState());
        }
        else if (class_name == "ClockState") {
            if (clock_defined) {
                std::cerr << "Clock State already defined for Agent '"
                            << sat->GetName() << "'\n";
                return;
            }
            CreateClockState(state_node, sat);
            clock_defined = true;
            Ptr<ClockState> clock_state_ptr = MakePtr<ClockState>(sat->GetClockState());

            InsertMap<IState>(state_map_, state_name, clock_state_ptr);
        }
        // ToDO: Add AttitudeState
        else {
            std::cerr << "Invalid State type: " << prior_field + ".class = " << class_name << "\n";
        }
    }

    //==============================================================================
    // LoadSatelliteConfig
    //==============================================================================

    void ConfigReader::LoadSatelliteConfig(YAML::Node satellite_node, bool print_val) {
        // First create the agents (Spacecraft objects)
        for (YAML::const_iterator it = satellite_node.begin(); it != satellite_node.end(); ++it) {
            std::string name = it->first.as<std::string>();
            Ptr<Spacecraft> sat = MakePtr<Spacecraft>();
            sat->SetName(name);
            satellites_.push_back(sat);
        }

        int idx = 0;
        for (YAML::const_iterator it = satellite_node.begin(); it != satellite_node.end(); ++it) {

            YAML::Node child_node = it->second;

            satellites_[idx]->SetEpoch(time_config_.epoch0);

            bool orbit_defined = false;
            bool clock_defined = false;

            // states
            if (child_node["states"]) {
                for (YAML::const_iterator sit = child_node["states"].begin(); 
                     sit != child_node["states"].end(); 
                     ++sit) 
                {
                    std::string state_name = sit->first.as<std::string>();
                    YAML::Node state_node = sit->second;
                    CreateSatState(state_node, state_name, satellites_[idx],
                                   orbit_defined, clock_defined);
                }
            }

            // device
            if (child_node["device"]) {
                for (YAML::const_iterator devIt = child_node["device"].begin();
                     devIt != child_node["device"].end(); 
                     ++devIt)
                {
                    std::string device_name = devIt->first.as<std::string>();
                    YAML::Node device_node  = devIt->second;
                    CreateSatCommDevice(device_node, device_name, satellites_[idx]);
                }
            }

            // Insert to map (Spacecraft)
            InsertMap<Spacecraft>(satellites_map_, satellites_[idx]->GetName(), satellites_[idx]);
            idx++;
        }

        // Print values if requested
        if (print_val) {
            std::cout << "Agent Config:\n";
            for (auto &sat : satellites_) {
                std::cout << "  Agent Name: " << sat->GetName() << "\n";
            }
        }
    }

} // namespace lupnt