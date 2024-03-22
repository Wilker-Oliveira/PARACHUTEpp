/** @brief Implementation of the Process Elementary Properties class. */

#ifndef PROCESSELEMETARYPROPERTIES_CUDA_HPP
#define PROCESSELEMETARYPROPERTIES_CUDA_HPP

#include <array>
#include <vector>
#include <numbers>
#include <cmath>

#include "PEP_kernels.hpp"

template <short N, typename fpt>
/**
 * @brief The ProcessElementaryProperties class
 * implements the elementary properties of the Sum of Sinusoids model, described
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 4, section 4.2.
 */
class ProcessElementaryProperties {
public:
  int ThperBl = 512; /**< Number of threads per block.*/

  /**
   * @brief The CalcMeanvalue() function computes the mean value of the process.
   * @return the value 0, theoretical mean of the gaussian random process.
   */
  fpt CalcMeanValue() const;

  /**
   * @brief The CalcMeanPower() function computes the mean power of the process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @return the theoretical mean power of the process.
   */
  fpt CalcMeanPower(const std::array<fpt, N> &pathGains) const;

  /**
   * @brief The CalcACF() function computes the autocorrelation of the process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
   * @param tau a fpt, the parameter of the autocorrelation.
   * @return the theoretical autocorrelation value for a determined tau.
   */
  fpt CalcACF(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies, fpt tau) const;

  /**
   * @brief The CalcPSD() function computes the Power Spectral Density of the
process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @param dopplerFrequencies a fpt array, a array containing the Doppler
frequencies of the process.
   * @param var_f a fpt, the parameter of the PSD (frequency domain).
   * @return the theoretical PSD value for a determined var_f.
   */
  fpt CalcPSD(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies, fpt var_f, fpt eps) const;

  /**
   * @brief The CalcDiscreteDopplerSpectrum() function computes the discrete
Doppler Spectrum of the process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
    * @return a vector containing the Doppler Spectrum values.
   */
  std::vector<fpt> CalcDiscreteDopplerSpectrum(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies) const;

  /**
   * @brief The CalcDopplerSpread() function computes the Doppler spread of the
process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
   * @return the Doppler spread of the process.
   */
  fpt CalcDopplerSpread(const std::array<fpt, N> &pathGains,
                        const std::array<fpt, N> &dopplerFrequencies) const;

  /**
   * @brief The CalcPeriodicity() function computes the periodicity of the
process.
   * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
   * @return if exists, returns the periodicity, if not, returns 0.
   */
  fpt CalcPeriodicity(const std::array<fpt, N> &dopplerFrequencies) const; // TBD

  /**
   * @brief the CalcPDF() function computes the Probability Density of the
process.
   * @param pathGains a fpt array, the array containing the path gains of the
process.
   * @param x a fpt, the targeted value of the process.
   * @return the PDF of the process.
   */  
  fpt CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/) const;

};

template <short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanValue() const {
    return 0;
}

template <short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanPower(const std::array<fpt,N> &pathGains) const{
    // Creating the data vector pointers.
    fpt *d_vec, *h_res, ret;

    // Allocating the device data vector.
    cudaMalloc(&d_vec, N*sizeof(fpt));

    // Copying data from host to device.
    cudaMemcpy(d_vec, pathGains.data(), N*sizeof(fpt), cudaMemcpyHostToDevice);

    // Calculating the Mean Power.
    CalcMeanPowerKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_vec);
    cudaDeviceSynchronize();

    summationKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_vec);
    cudaDeviceSynchronize();    

    // Copying result from device to host.
    cudaMemcpy(h_res, d_vec, sizeof(fpt), cudaMemcpyDeviceToHost);

    // Copying result from h_res to ret.    
    ret = *h_res; 
    
    // Freeing used memory.
    cudaFree(d_vec);
    delete h_res;

    return ret;
}

template <short N, typename fpt>
fpt ProcessElementaryProperties<N, fpt>::CalcACF(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies, fpt tau) const {
    // Creating the data vector pointers.
    fpt *d_pG, *d_dF, *d_ACF, *h_res, ret;

    // Allocating the device data vector.
    cudaMalloc(&d_pG, N * sizeof(fpt));
    cudaMalloc(&d_dF, N * sizeof(fpt));
    cudaMalloc(&d_ACF, N * sizeof(fpt));

    // Copying data from host to device.
    cudaMemcpy(d_pG, pathGains.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dF, dopplerFrequencies.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);

    // Calculating the ACF.
    CalcACFKernel<N,fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_pG, d_dF, d_ACF, tau);
    cudaDeviceSynchronize();

    summationKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_ACF);
    cudaDeviceSynchronize();

    // Copying result from device to host.
    cudaMemcpy(h_res, d_ACF, sizeof(fpt), cudaMemcpyDeviceToHost);

    // Copying result from h_res to ret.    
    ret = *h_res;

    // Freeing used memory.
    cudaFree(d_pG);
    cudaFree(d_dF);
    cudaFree(d_ACF);
    delete h_res;

    return ret;
}

template <short N, typename fpt>
fpt ProcessElementaryProperties<N, fpt>::CalcPSD(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies, fpt var_f, fpt eps = 0.0001) const {
    // Creating the data vector pointers.
    fpt *d_pG, *d_dF, *d_Sf, *h_res, ret;

    // Allocating the device data vector.    
    cudaMalloc(&d_pG, N * sizeof(fpt));
    cudaMalloc(&d_dF, N * sizeof(fpt));
    cudaMalloc(&d_Sf, N * sizeof(fpt));

    // Copying data from host to device.    
    cudaMemcpy(d_pG, pathGains.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dF, dopplerFrequencies.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);

    // Calculating the PSD.
    CalcPSDKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_pG, d_dF, d_Sf, var_f, eps);
    cudaDeviceSynchronize();

    summationKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_Sf);
    cudaDeviceSynchronize();

    // Copying result from device to host.    
    cudaMemcpy(h_res, d_Sf, sizeof(fpt), cudaMemcpyDeviceToHost);

    // Copying result from h_res to ret.    
    ret = *h_res;

    // Freeing used memory.    
    cudaFree(d_pG);
    cudaFree(d_dF);
    cudaFree(d_Sf);
    delete h_res;

    return ret;
}

template <short N, typename fpt>
std::vector<fpt>
ProcessElementaryProperties<N, fpt>::CalcDiscreteDopplerSpectrum(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies) const {
    // Creating the data vector pointers.
    fpt *d_pG, *d_dF, *d_dDS, *h_res;
    std::vector<fpt> ret(4*N,0);

    // Allocating the device data vector.    
    cudaMalloc(&d_pG, N * sizeof(fpt));    
    cudaMalloc(&d_dF, N * sizeof(fpt));
    cudaMalloc(&d_dDS, 4 * N * sizeof(fpt));

    // Allocating the host data vector.    
    h_res = new fpt[4*N]();

    // Copying data from host to device.    
    cudaMemcpy(d_pG, pathGains.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dF, dopplerFrequencies.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);

    // Calculating the discrete Doppler spectrum.    
    CalcDDSKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_pG, d_dF, d_dDS);
    cudaDeviceSynchronize();

    // Copying result from device to host.    
    cudaMemcpy(h_res, d_dDS, 4 * N * sizeof(fpt), cudaMemcpyDeviceToHost);

    // Copying result from h_res to ret.    
    for (int i = 0; i < 4 * N; i++) {
      ret[i] = h_res[i];
    }

    // Freeing used memory.    
    cudaFree(d_pG);
    cudaFree(d_dF);
    cudaFree(d_dDS);
    delete[] h_res;

    return ret;
}

template <short N, typename fpt>
fpt ProcessElementaryProperties<N, fpt>::CalcDopplerSpread(const std::array<fpt, N> &pathGains, const std::array<fpt, N> &dopplerFrequencies) const {
    // Creating the data vector pointers.
    fpt *d_pG, *d_dF, *d_DS, *h_res, ret, sigma=std::sqrt(CalcMeanPower(pathGains));
    fpt sum, beta_i;

    // Allocating the device data vector.    
    cudaMalloc(&d_pG, N * sizeof(fpt));
    cudaMalloc(&d_dF, N * sizeof(fpt));
    cudaMalloc(&d_DS, N * sizeof(fpt));

    // Copying data from host to device.    
    cudaMemcpy(d_pG, pathGains.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dF, dopplerFrequencies.data(), N * sizeof(fpt), cudaMemcpyHostToDevice);

    // Calculating the Doppler Spread.
    CalcDopplerSpreadKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_pG, d_dF, d_DS);
    cudaDeviceSynchronize();

    summationKernel<N, fpt><<<(N + ThperBl - 1) / ThperBl, ThperBl>>>(d_DS);
    cudaDeviceSynchronize();

    // Copying data from device to host.    
    cudaMemcpy(h_res, d_DS, sizeof(fpt), cudaMemcpyDeviceToHost);

    // Storing results in ret.    
    sum = *h_res;
    beta_i = 2*pow(M_PI,2)*sum;
    ret = sqrt(beta_i)/(2*M_PI*sigma);

    // Freeing used memory.    
    cudaFree(d_pG);
    cudaFree(d_dF);
    cudaFree(d_DS);
    delete h_res;
    return ret;
}

template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPeriodicity(const std::array<fpt,N> &dopplerFrequencies) const {
}// finish later

template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/) const {

  if(fabs(x) < pathGain){
    return 1/(M_PI*fabs(pathGain)*pow(1-pow(x/pathGain,2),1/2) );
  }
  else return 0;
}
   


#endif
