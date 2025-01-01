/**
 * @file GnssTransmitter.cpp
 * @author Stanford NAV LAB
 * @brief Handles Gnss transmit
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/gnss_transmitter.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/file.h"
#include "lupnt/data/kernels.h"
#include "lupnt/measurements/gnss_channel.h"
#include "lupnt/numerics/string_utils.h"

namespace lupnt {

  bool StringInCandidates(std::vector<std::string>& candidates, std::string str) {
    for (size_t i = 0; i < candidates.size(); i++) {
      if (str.find(candidates[i]) != std::string::npos) {
        return true;
      }
    }
    return false;
  }

  void GnssTransmitter::InitializeGnssTransmitter() {
    // std::cout << "GNSS Transmitter: " << gnss_type_ << std::endl;

    if (gnss_type_ == GnssType::GPS) {
      InitializeGPSTransmitter();
      gnss_type_str_ = "GPS";
    } else if (gnss_type_ == GnssType::GLONASS) {
      InitializeGLONASSTransmitter();
      gnss_type_str_ = "GLONASS";
    } else if (gnss_type_ == GnssType::GALILEO) {
      InitializeGALILEOTransmitter();
      gnss_type_str_ = "GALILEO";
    } else if (gnss_type_ == GnssType::BEIDOU) {
      InitializeBEIDOUTransmitter();
      gnss_type_str_ = "BEIDOU";
    } else if (gnss_type_ == GnssType::QZSS) {
      InitializeQZSSTransmitter();
      gnss_type_str_ = "QZSS";
    } else {
      throw std::runtime_error("Invalid GNSS type");
    }
  }

  /**
   * @brief Initialize the GPS transmitter
   *
   */
  void GnssTransmitter::InitializeGPSTransmitter() {
    // Todo: Autonomously generate the prn to svn mapping
    std::filesystem::path csvpath(GetDataPath() / "gnss" / "gps_table.csv");
    std::vector<std::vector<std::string>> gps_table = ReadCSV(csvpath.string());

    for (size_t i = 0; i < gps_table.size(); i++) {
      if (std::stoi(gps_table[i][0]) == prn_) {
        std::string gps_type = gps_table[i][2];  // 'IIA', 'IIR', 'IIR-M', 'IIF'

        // set transmittion power and antenna pattern, depending on the gps type
        std::string ant_name;
        if (gps_type == "IIA") {
          P_tx = 14.3;                 // dB-W
          ant_name = gps_table[i][4];  // ACE Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIR") {
          P_tx = 15.0;                 // dB_W
          ant_name = gps_table[i][3];  // LM Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIR-M") {
          P_tx = 15.0;                 // dB_W
          ant_name = gps_table[i][3];  // LM Pattern
          freq_list = {"L1", "L2"};
        } else if (gps_type == "IIF") {
          P_tx = 14.3;                 // dB_W
          ant_name = gps_table[i][4];  // ACE Pattern
          freq_list = {"L1", "L2", "L5"};
        } else if (gps_type == "III") {
          P_tx = 14.3;                 // dB_W
          ant_name = gps_table[i][3];  // LM Pattern
          freq_list = {"L1", "L2", "L5"};
        } else {
          std::runtime_error("Invalid GPS type");
        }

        // std::cout << "PRN: " << prn_ << " type: " << gps_type << " Antenna: " << ant_name
        //           << std::endl;

        if (gps_type == "III"){   // For Block III, the antenna pattern is different for each frequency
          for (auto freq : freq_list) {
            antenna_[freq] = Antenna(ant_name + "_" + freq + ".txt");
          }
        }
        else {  // For other GPS types, we assume the antenna pattern is the same for all frequencies
          for (auto freq : freq_list) {
            antenna_[freq] = Antenna(ant_name);
          }
        }

        break;
      }
    }
    return;
  }

  /**
   * @brief   Initialize the GLONASS transmitter
   *
   */
  void GnssTransmitter::InitializeGLONASSTransmitter() {
    std::cout << "Antenna type not implemented yet for " << gnss_type_str_ << std::endl;
  }

  /**
   * @brief Initialize the GALILEO transmitter
   *
   */
  void GnssTransmitter::InitializeGALILEOTransmitter() {
    freq_list = {"E1", "E5a", "E5b", "E6"};
    for (auto freq : freq_list) {
      antenna_[freq] = Antenna("Galileo_" + freq + ".txt"); 
    }
    
    // The antenna pattern given from ESA is the EIRP=gain + P_tx
    // Since we do not know the trans power, we assumed it to be 14.0 dBW when we generated the data
    // Therefore we set the P_tx to 14.0 dBW
    P_tx = 14.0;   // dB-W  Assume fixed
  }

  /**
   * @brief Initialize the BEIDOU transmitter
   */
  void GnssTransmitter::InitializeBEIDOUTransmitter() {
    std::cout << "Antenna type not implemented yet for " << gnss_type_str_ << std::endl;
  }

  /**
   * @brief 
   * 
   */
  void GnssTransmitter::InitializeQZSSTransmitter() {
    std::vector<std::string> names = {"1R", "02", "03", "04", "05", "06", "07"};

    if (prn_ <= 4) {
      freq_list = {"L1", "L2", "L5"};
    } else {
      freq_list = {"L1", "L5"};
    }

    for (auto freq : freq_list) {
      antenna_[freq] = Antenna("QZSS_" + names[prn_-1] + "_" + freq + ".txt");
    }

    // Reference: Enhancing Navigation Accuracy in a Geostationary Orbit by
    //             Utilizing a Regional Navigation Satellite System
    P_tx = 14.1;  // dB-W

    return;

  }

  /**
   * @brief Get the transmitter orientation
   *        Todo: Change this depending on the gnss type
   *
   * @param t
   * @param r_tx_gcrf
   * @return Vec3d
   */
  std::vector<Vec3d> GnssTransmitter::GetTransmitterOrientation(double t, Vec3d& r_tx_gcrf) {
    // (Sun-Earth) - (Sat-Earth)
    Vec3d r_sat2sun
        = GetBodyPosVel(t, NaifId::EARTH, NaifId::SUN, Frame::GCRF).cast<double>().head(3)
          - r_tx_gcrf;
    auto e_z_gnss = -r_tx_gcrf.normalized();  // Face towards earth center
    auto e_y_gnss = r_sat2sun.cross(r_tx_gcrf).normalized();
    auto e_x_gnss = e_y_gnss.cross(e_z_gnss).normalized();

    std::vector<Vec3d> e_gnss = {e_x_gnss, e_y_gnss, e_z_gnss};
    return e_gnss;
  }


  double GnssTransmitter::GetTransmitterAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) {
    // Get the first freq in antenna map
    std::string freq = freq_list[0];
    double At = GnssTransmitter::GetTransmitterAntennaGainFreq(t, r_tx_gcrf, r_rx_gcrf, freq); // use L1 for default
    return At;
  }

  double GnssTransmitter::GetTransmitterAntennaGainFreq(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf, std::string freq) {
    auto e_gnss = GnssTransmitter::GetTransmitterOrientation(t, r_tx_gcrf);
    Vec3d e_x_gnss = e_gnss[0];
    Vec3d e_y_gnss = e_gnss[1];
    Vec3d e_z_gnss = e_gnss[2];
    Vec3d u_tx_rx = (r_rx_gcrf - r_tx_gcrf).normalized();
    double phi_tx = acos(u_tx_rx.dot(e_z_gnss));
    double theta_tx = atan2(u_tx_rx.dot(e_y_gnss), u_tx_rx.dot(e_x_gnss));
    double At = GnssTransmitter::ComputeGain(theta_tx, phi_tx, freq).val();   // use L1 for default
    return At;
  }

  /**
   * @brief Generate a transmission
   *
   * @param t
   * @return Transmission
   */
  GnssTransmission GnssTransmitter::GenerateTransmission(double t) {
    CartesianOrbitState cart_state = GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    ConvertOrbitStateFrame(cart_state, t, Frame::GCRF);

    GnssTransmission trans;
    trans.dt_tx = 0.0;
    trans.r_tx = cart_state.r().cast<double>();
    trans.v_tx = cart_state.v().cast<double>();
    return trans;
  }

}  // namespace lupnt
