/**
 * @file comm_utils.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-08-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once
#include "lupnt/core/definitions.h"

namespace lupnt {

  enum Modulation {
    Residual,  // Residual carrier
    // Suppressed carrier
    BPSK,     // Binary Phase Shift Keying
    QPSK,     // Quadrature Phase Shift Keying
    OQPSK,    // Offset Quadrature Phase Shift Keying
    GMSK,     // Gaussian Minimum Shift Keying
    GMSK_PN,  // GMSK with PN modulation
  };

  enum FrequencyBand { UHF, L, S, Cband, X, Ku, K, Ka };

  /**
   * @brief Compute the Approximate BER for different Carrier modulation
   *  https://www.unilim.fr/pages_perso/vahid/notes/ber_awgn.pdf
   * @param PT_N0
   * @param Modulation
   * @return Real
   */
  Real ComputeBER(Real EbN0, Modulation modulation_type);

  /**
   * @brief Compute Es/N0 from Eb/N0
   *
   * @param EbN0  Energy per bit to noise spectral density ratio
   * @param modulation_order  Modulation order (QPSK: 2, BPSK: 1)
   * @param coding_rate  Coding rate
   *
   */
  Real ComputeEsN0(Real EbN0, Real modulation_order, Real coding_rate);

  /**
   * @brief Compute the number of bits per symbol
   *
   * @param modulation_type  Modulation type
   * @return Real
   */
  Real BitsPerSymbol(Modulation modulation_type);
  /**
   * @brief Get the Frequency Band object
   *
   * @param f   frequency [Hz]
   * @return FrequencyBand
   */
  FrequencyBand GetFrequencyBand(Real f);

  /**
   * @brief Get the (recommended) Transponder Turn Around Ratio for spacecraft
   *  // https://deepspace.jpl.nasa.gov/dsndocs/810-005/201/201B.pdf
   *
   * @param fbu  uplink frequency band
   * @param fbd  downlink frequency band
   * @return Real   Turn around ratio
   */
  Real GetTransponderTurnAroundRatio(FrequencyBand fbu, FrequencyBand fbd);

  /**
   * @brief Compute the carrier loop signal-to-noise ratio
   *
   * @param PT_N0  downlink total signal power to noise spectral density ratio
   * [Hz]
   * @param B_L_carrier  carrier loop noise bandwidth [Hz]
   * @param T_s           period of the binary symbol [s]
   * @param modulation_type  modulation type
   * @param m_R           modulation index (only for residual carrier)
   *
   */
  Real ComputeCarrierLoopSNR(Real PT_N0, Real B_L_carrier, Real T_s, Modulation modulation_type,
                             Real m_R);

}  // namespace lupnt
