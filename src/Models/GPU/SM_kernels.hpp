#ifndef SM_KERNELS_H
#define SM_KERNELS_H

//Cuda libraries
#include <cuda.h>
#include <cuda_runtime.h>

template <short N, typename fpt>
/**
 * @brief The CalcProcessKernel() calculates the response of the Gaussian
process in time.
 * @param time a fpt vector, the time vector.
 * @param pathGains a fpt array, the array containing the path gains of the
process.
 * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
 * @param processSample a fpt array, the output data vector pointer.
 * @param NofSamples a fpt, the number of samples to be computed.
 */
__global__ void  CalcProcessKernel(fpt* time, fpt* pathGains, fpt* dopplerFrequencies, fpt* phases, fpt* processSample, int NofSamples) {

  int j = threadIdx.x + blockDim.x * blockIdx.x;

  if(j<NofSamples){
    for(int i=0;i<N;i++)
      processSample[j]+=pathGains[i]*cos(2*M_PI*dopplerFrequencies[i]*time[j] + phases[i]);
  }
}

template <short N, typename fpt>
/**
 * @brief The CalcProcessKernel() calculates the response of the Hilbert
process in time.
 * @param time a fpt vector, the time vector.
 * @param pathGains a fpt array, the array containing the path gains of the
process.
 * @param dopplerFrequencies a fpt array, the array containing the Doppler
frequencies of the process.
 * @param processSample a fpt array, the output data vector pointer.
 * @param NofSamples a fpt, the number of samples to be computed.
 */

__global__ void  CalcHilbertProcessKernel(fpt* time, fpt* pathGains, fpt* dopplerFrequencies, fpt* phases, fpt* processSample, int NofSamples) {

  int j = threadIdx.x + blockDim.x * blockIdx.x;

  if(j<NofSamples){
    for(int i=0;i<N;i++)
      processSample[j]+=pathGains[i]*sin(2*M_PI*dopplerFrequencies[i]*time[j] + phases[i]);
  }

}

#endif
