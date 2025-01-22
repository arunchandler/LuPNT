/**
 * @file gnss_transmitter.h
 * @author Stanford NAV LAB
 * @brief  Handles Gnss transmit
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string>
#include <tuple>

#include "antenna.h"
#include "gnss_measurement.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/transmission.h"

namespace lupnt {

  enum GnssType {
    GPS,      // Global Positioning System
    GLONASS,  // Global Navigation Satellite System
    GALILEO,  // European Global Navigation Satellite System
    BEIDOU,   // Chinese Global Navigation Satellite System
    QZSS,     // Quasi-Zenith Satellite System
  };

  class GnssChannel;

  class GnssTransmitter : public Transmitter {
  public:
    std::map<std::string, Antenna> antenna_;  // Antenna gain pattern [deg & dB]
    std::string gnss_type_str_;               // Name of the atenna system
    GnssType gnss_type_;                      // Type of the GNSS system
    std::string txrx = "TX";                  // Type of comms system
    int prn_;                                 // PRN of the transmitter satellite
    double Rc;                                // Ranging chip rate [Hz]
    std::vector<std::string> freq_list;       // List of frequencies (by signal names)
    std::map<std::string, double> freq_map;   // map string to frequencies
    std::map<std::string, double> rc_map;     // map string to chip rates

    // Tramsmitter
    GnssTransmitter(GnssType gnss_type, int prn) : gnss_type_(gnss_type), prn_(prn) {
      freq_map = {{"L1", 1575.42e6},  {"L2", 1227.60e6}, {"L5", 1176.45e6},
                  {"E1", 1575.42e6},  {"E6", 1278.75e6}, {"E5", 1191.795e6},
                  {"E5a", 1176.45e6}, {"E5b", 1207.14e6}};  // frequency maps

      // Information from https://gnss-sdr.org/docs/tutorials/gnss-signals/
      rc_map
          = {{"L1", 1.023e6}, {"L2", 0.5115e6}, {"L5", 10.23e6}, {"E1", 1.023e6}, {"E6", 0.5115e6},
             {"E5", 10.23e6}, {"E5a", 10.23e6}, {"E5b", 10.23e6}};  // chip rate maps

      InitializeGnssTransmitter();
    }

    // Transmitter Initializations
    void InitializeGnssTransmitter();
    void InitializeGPSTransmitter();
    void InitializeGLONASSTransmitter();
    void InitializeGALILEOTransmitter();
    void InitializeBEIDOUTransmitter();
    void InitializeQZSSTransmitter();

    // Get transmitter orientatiion
    std::vector<Vec3> GetTransmitterOrientation(Real t, const Vec3& rv_tx_gcrf);

    Real GetTransmitterAntennaGain(Real t, const Vec3& r_tx_gcrf, const Vec3& r_rx_gcrf) override;
    Real GetTransmitterAntennaGainFreq(Real t, const Vec3& r_tx_gcrf, const Vec3& r_rx_gcrf,
                                       const std::string& freq);

    // Get the Transmittion Information
    GnssTransmission GenerateTransmission(Real t);

    // Getters and Setters
    void SetChannel(const Ptr<SpaceChannel> ch) override {
      channel_ = std::static_pointer_cast<GnssChannel>(ch);
    };
    int GetPRN() { return prn_; };
    GnssType GetGnssType() { return gnss_type_; };
    Real ComputeGain(Vec3d direction, const std::string& freq) {
      return antenna_[freq].ComputeGain(direction);
    };
    Real ComputeGain(Real theta, Real phi, const std::string& freq) {
      return antenna_[freq].ComputeGain(theta, phi);
    };

  private:
    Ptr<GnssChannel> channel_;  // Channel that the device is connected to
  };
}  // namespace lupnt
