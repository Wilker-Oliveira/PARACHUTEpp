#ifndef PEP_KERNELS_H
#define PEP_KERNELS_H

//Cuda libraries
#include <cuda.h>
#include <cuda_runtime.h>
//Kernels
template <short N, typename fpt>
/**
 * @brief The kernel summationkernel() performs a parallel sum of the elements
in a 1-D array.
 * @param d_input a fpt, the data vector pointer with the elements to be summed.
 */
__global__ void summationKernel(fpt *d_input) {

  const int id = threadIdx.x + blockDim.x * blockIdx.x;
 
  auto step = 1;
  int num_threads = blockDim.x * gridDim.x;

  while (num_threads > 0) {

    if (id < num_threads) {
      const auto fst = id * step * 2;
      const auto scd = fst + step;

      d_input[fst] += d_input[scd];
    }
    step <<= 1;
    num_threads >>=1;    
  }   
}

template <short N, typename fpt>
/**
 * @brief The CalcMeanPowerKernel() computes the elements to be added in order
to calculate the theoretical mean power of the process, shown in equation (4.10) in "Mobile Radio Channels", by Matthias Patzold.
 * @param d_input a fpt, the data vector pointer with the path gains of the process.
 */
__global__ void CalcMeanPowerKernel(fpt *d_input) {

  int id = threadIdx.x + blockDim.x * blockIdx.x;
  if (id < N) {
    d_input[id] = pow(d_input[id],2)/2;
  }
}

template <short N, typename fpt>
/**
 * @brief The CalcACFKernel() computes the elements to be added in the summation of the autocorrelation expression in equation (4.11) in "Mobile Radio Channels", by Matthias Patzold.
 * @param d_pg a fpt, the data vector pointer with the path gains of the
process.
 * @param d_df a fpt, the data vector pointer with the Doppler frequencies of
the process.
 * @param d_out a fpt, the output data vector pointer.
 * @param tau a fpt, the parameter of the autocorrelation.
 */
__global__ void CalcACFKernel(fpt *d_pg, fpt *d_df, fpt *d_out, fpt tau) {

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < N) {
      d_out[id] = pow(d_pg[id],2)*cos(2*M_PI*d_df[id]*tau)/2;
    }
  }

template <short N, typename fpt>
/**
 * @brief The CalcPSDKernel() computes the elements to be added in the summation of power spectral density expression shown in equation (4.18) in "Mobile Radio Channels", by Matthias Patzold, for a single frequency value.
 * @param d_pg a fpt, the data vector pointer with the path gains of the
process.
 * @param d_df a fpt, the data vector pointer with the Doppler frequencies of
the process.
 * @param d_out a fpt, the output data vector pointer.
 * @param var a fpt, the parameter of the PSD (frequency domain).
 * @param eps a fpt, the maximum allowed gap between the Doppler frequency and the var_f values.
 */  
__global__ void CalcPSDKernel(fpt * d_pg, fpt * d_df, fpt *d_out, fpt var, fpt eps) {

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < N) {
      if (fabs(var - d_df[id]) <= eps || fabs(d_df[id] + var) <= eps) {
        d_out[id] = pow(d_pg[id],2)/4;
      } else {
        d_out[id] = 0;
      }     
    }
  }

template <short N, typename fpt>
/**
 * @brief The CalcDDSKernel() computes the discrete Doppler spectrum of the process.
 * @param d_pg a fpt, the data vector pointer with the path gains of the
process.
 * @param d_df a fpt, the data vector with the Doppler frequencies of the
process.
 * @param d_out a fpt, the output data vector pointer.
 */  
__global__ void CalcDDSKernel(fpt * d_pg, fpt * d_df, fpt * d_out) {

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < N) {
      d_out[2 * id] = -d_df[id];
      d_out[(2 * id) + 1] = pow(d_pg[id], 2) / 4;
      d_out[(2 * id) + (2 * N)] = d_df[id];
      d_out[(2 * id) + (2 * N) + 1] = pow(d_pg[id], 2) / 4;      
    }
  }

template <short N, typename fpt>
/**
* @brief The CalcDopplerSpreadKernel() computes the elements in the summation of
the expression of the Doppler spread of the process, in equation (4.29) in
"Mobile Radio Channels", by Matthias Patzold.
 * @param d_pg a fpt, the data vector pointer with the path gains of the
process.
 * @param d_df a fpt, the data vector pointer with the Doppler frequencies of
the process.
 * @param d_out a fpt, the output data vector pointer.
 */  
__global__ void CalcDopplerSpreadKernel(fpt * d_pg, fpt * d_df, fpt * d_dout) {

    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < N) {
      d_dout[id] = pow(d_pg[id]*d_df[id],2);
    }
  }  



#endif
