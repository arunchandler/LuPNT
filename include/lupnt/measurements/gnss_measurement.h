/**
 * @file gnss_measurement.h
 * @author Stanford NAV LAB
 * @brief Class that constructs GPS measurement data from channel observables
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <memory>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/measurements/transmission.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  enum GnssMeasurementType {
    PR,   // Pseudorange
    PRR,  // Pseudorange rate
    CP,   // Carrier phase
  };

  class GnssMeasurement {
    // Implemenation based on Gnss SDR Observables block:
    // https://gnss-sdr.org/docs/sp-blocks/observables/
  private:
    int n_meas = 0;                             // Number of measurements
    std::vector<GnssTransmission> trans_store;  // list of transmittion data
    std::vector<int> ID_tx;                     // ID of the transmitter (n_meas)
    std::vector<std::string> gnss_types_;       // GNSS type (n_meas)

    // Visbility
    VecXd vis_earth;
    VecXd vis_moon;
    VecXd vis_antenna;
    VecXd vis_atmos;
    VecXd vis_ionos;

    // Pseudorange measurement
    Real t_rx;  // Signal reception time measured by the receiver clock [s]
    Real t_tx;  // Signal transmission time measured by the transmitter clock [s]

    VecX P_rx;    // Pseudorange measurement [km] (n_meas * n_bands)
    VecX rho_rx;  // True range from the transmitter to the receiver’s
    // antenna [km] (n_meas)

    Real dt_rx;  // Receiver clock offset from Gnss time [s]
    VecX dt_tx;  // Transmitter clock offset from Gnss time [s] (n_meas)

    VecX I_rx;   // Ionospheric delay [km] (n_meas * n_bands)
    VecX T_rx;   // Tropospheric delay [km] (n_meas)
    VecX eps_P;  // Pseudorange measurement noise [km] (n_meas)

    // Carrier phase measurement
    VecX phi_rx;     // Carrier phase measurement [cycles] (n_meas * n_bands)
    VecX phi_rx_tx;  // Phase of the receiver's local oscillator at
                     // time t [cycles] (n_meas * n_bands)
    VecX phi_tx;     // Phase of the transmitted signal at time t [cycles]
                     // (n_meas * n_bands)

    VecX N_rx;      // Carrier phase int ambiguity [cycles] (n_meas * n_bands)
    VecX eps_phi;   // Carrier phase measurement noise [cycles]
                    // (n_meas * n_bands)
    VecX f;         // Carrier frequency [Hz] (n_bands) // TODO: Convert to vector
    VecX phi_0;     // Initial phase of the transmitted signal at time t_0
                    // [cycles] (n_meas * n_bands)
    VecX lambda_;   // Carrier wavelength [km] (n_bands) // TODO: Convert to vector
    VecX L_rx_txr;  // Integer component of the receiver’s numerically
                    // controlled oscillator (NCO) phase
    VecX K_rx_txr;  // Integer component of the propagation term

    // Phase-range masurement
    VecX Phi_rx;     // Phase-range measurement [km] (n_meas * n_bands)
                     // (n_meas * n_bands)
    VecX B_rx;       // Carrier phase bias [cycles] (n_meas * n_bands)
    VecX dPhi_rx;    // Receiver’s antenna phase center variation [km]
                     // (n_meas * n_bands)
    VecX d_rx_pco;   // Receiver’s antenna phase center offset in local
                     // coordinates [km] (n_meas) (n_bands x 3)
    VecX d_rx_pcv;   // Receiver’s antenna phase center variation [km]
                     // (n_meas * n_bands)
    VecX d_tx_pco;   // Transmitter’s antenna phase center offset in
                     // local coordinates [km] (n_meas * n_bands x 3)
    VecX d_tx_pcv;   // Transmitter’s antenna phase center variation
                     // [km] (n_meas * n_bands)
    MatX e_rx_enu;   // LOS vector from receiver antenna to satellite in
                     // local coordinates [km] (n_meas x 3)
    MatX e_rx;       // LOS vector from receiver antenna to satellite in
                     // ECEF coordinates [km] (n_meas x 3)
    VecX E;          // Coordinates transformation matrix from the satellite
                     // body‐fixed coordinates to ECEF coordinates
    Vec3 d_rx_disp;  // Displacement by Earth tides at the receiver
                     // position in local coordinates [km] (3)
    VecX eps_Phi;    // Phase-range measurement noise [km] (n_bands)

    // Doppler shift measurement
    VecX D_rx;   // Doppler shift measurement [Hz] (n_meas * n_bands)
    VecX r_rx;   // Position of the receiver at time t_rx [km] (3)
    MatX r_tx;   // Position of the transmitter at time t_rx [km] (n_meas x 3)
    VecX v_rx;   // Velocity of the receiver at time t_rx [m/s] (3)
    MatX v_tx;   // Velocity of the transmitter at time t_tx [m/s] (n_meas x 3)
    VecX eps_D;  // Doppler shift measurement noise [Hz] (n_meas * n_bands)

    // Pseudorange rate measurement
    VecX PR_rx;      // Pseudorange rate measurement [m/s] (n_meas * n_bands)
    Real dt_rx_dot;  // Receiver clock drift [s/s]
    VecX dt_tx_dot;  // Transmitter clock drift [s/s] (n_meas)
    VecX eps_PR;     // Pseudorange rate measurement noise [m/s]
                     // (n_meas * n_bands)

    // Receiver param
    GnssReceiverParam gnssr_param;
    Real chip_rate;  // Chip rate [Hz]

    // Link budget
    VecX CN0;  // Carrier‐to‐noise density [dB‐Hz] (n_meas * n_bands)

  public:
    GnssMeasurement(const std::vector<GnssTransmission> transmissions);

    GnssMeasurement ExtractSignal(std::string freq_label);
    GnssMeasurement ExtractSignal(std::vector<std::string> freq_labels);

    GnssMeasurement ApplyIonoMask();

    // Transmission data
    int GetTrackedSignalNum() const { return n_meas; }

    std::vector<int> GetTxIds() const { return ID_tx; }
    VecX GetCN0() const { return CN0; }
    VecXd GetEarthOccultation() const { return vis_earth; }
    VecXd GetMoonOccultation() const { return vis_moon; }
    VecXd GetAntennaOccultation() const { return vis_antenna; }
    VecXd GetIonosOccultation() const { return vis_ionos; }
    VecXd GetAtmosOccultation() const { return vis_atmos; }

    /***********************************************************
     * General Methods for computing Measurements
     ***********************************************************/

    /**
     * @brief Compute the pseudorange measurement
     *
     * @param r_rx   Receiver position [km]
     * @param dt_rx  Receiver clock offset from Gnss time [s]
     * @param with_noise   use noise
     * @param seed   random seed
     * @return VecX
     */
    VecX ComputeGnssPseudorange(const VecX& r_rx, Real dt_rx, bool with_noise = false,
                                int seed = 0);
    VecX ComputeGnssPseudorangerate(const VecX& r_rx, const VecX& v_rx, Real dt_rx_dot,
                                    bool with_noise = false, int seed = 0);
    VecX ComputeGnssCarrierPhase(const VecX& r_rx, Real dt_rx, const VecX& N_rx,
                                 bool with_noise = false, int seed = 0);

    /***********************************************************
     *  Methods for true measurement generation
     ***********************************************************/

    /**
     * @brief Get the Gnss Measurement for the observed signal
     *
     * @param meas_type  vector of measurement types
     * @param with_noise  use noise
     * @param seed   random seed
     * @return VecX
     */
    VecX GetGnssMeasurement(std::vector<GnssMeasurementType> meas_type, bool with_noise = false,
                            int seed = 0);
    /**
     * @brief Get the Pseudorange for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetPseudorange(bool with_noise = false, int seed = 0);

    /**
     * @brief Get the Pseudorange rate for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetPseudorangerate(bool with_noise = false, int seed = 0);

    /**
     * @brief Get the Carrier Phase for the observed signal
     *
     * @param with_noise  use noise
     * @param seed  random seed
     * @return VecX
     */
    VecX GetCarrierPhase(bool with_noise = false, int seed = 0);

    /***********************************************************
     *  Methods for predicted measurement generation
     ***********************************************************/
    /**
     * @brief Get the Gnss Measurement for the predicted state
     *       This method is used for the measurement update step
     *       in the filter
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param frame_in  coordinate system of the input state
     */
    VecX GetPredictedGnssMeasurement(Real epoch, const Vec6& rv_pred, const Vec2& clk_pred,
                                     const VecX& N_pred, std::vector<GnssMeasurementType> meas_type,
                                     Frame frame_in = Frame::MOON_CI, MatXd* H_gnss = nullptr);

    /**
     * @brief Get the Pseudorange for the predicted state
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_pr       Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedPseudorange(Real epoch, const Vec6& rv_pred, const Vec2& clk_pred,
                                 Frame frame_in = Frame::MOON_CI, MatXd* H_pr = nullptr);

    /**
     * @brief Get the Pseudorange Analytical Jacobian object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_pr      Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return * VecX
     */
    VecX GetPredictedPseudorangeAnalyticalJacobian(Real epoch, const Vec6& rv_pred,
                                                   const Vec2& clk_pred,
                                                   Frame frame_in = Frame::MOON_CI,
                                                   MatXd* H_pr = nullptr);

    /**
     * @brief Get the Pseudorange Rate object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_prr     Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedPseudorangerate(Real epoch, const Vec6& rv_pred, const Vec2& clk_pred,
                                     MatXd& H_prr, Frame frame_in = Frame::MOON_CI);

    /**
     * @brief Get the Carrier Phase object
     *
     * @param epoch     epoch time
     * @param rv_pred   predicted position and velocity
     * @param clk_pred  predicted clock offset and drift
     * @param H_cp      Jacobian of the measurement function
     * @param frame_in  coordinate system of the input state
     * @return VecX
     */
    VecX GetPredictedCarrierPhase(Real epoch, const Vec6& rv_pred, const Vec2& clk_pred,
                                  VecX N_pred, MatXd& H_cp, Frame frame_in = Frame::MOON_CI);

    /*********************************************************************
     * Noise Models
     ********************************************************************/
    VecX GetGnssNoiseStdVec(std::vector<GnssMeasurementType> meas_type);
    VecX GetPseudorangeNoiseStdVec();
    VecX GetPseudorangeRateNoiseStdVec();
    VecX GetCarrierPhaseNoiseStdVec();

    void SetGnssReceiverParam(GnssReceiverParam gnssr_param_input) {
      gnssr_param = gnssr_param_input;
    }

    /**
     * @brief Compute the pseudorange noise using thermal noise in DLL
     * Reference: Reference: "Understanding GPS", p195
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return Real  Pseudorange noise [km]
     */
    Real ComputeGnssPseudorangeNoise(Real CN0);

    /**
     * @brief Compute the pseudorange rate noise using thermal noise in FLL
     * Reference: "Understanding GPS", p192
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return Real  Pseudorange rate noise [km/s]
     */
    Real ComputeGnssPseudorangerateNoise(Real CN0, Real lambda);

    /**
     * @brief Compute the carrier phase noise using thermal noise in PLL
     * Reference: "Understanding GPS", p185
     *
     * @param CN0  Carrier‐to‐noise density [dB‐Hz]
     * @return Real  Carrier phase noise [cycles]
     */
    Real ComputeGnssCarrierPhaseNoise(Real CN0, Real lambda);
  };
}  // namespace lupnt
