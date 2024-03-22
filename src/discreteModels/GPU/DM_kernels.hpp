#ifndef DM_KERNELS_H
#define DM_KERNELS_H

//Kernels
template <typename fpt>
/**
 * @brief The kernel rayleighDistributionKernel() computes samples of the Rayleigh distribution using the Percentile Transformation Method presented in "Probility, Random Variables and Stochastic Processes" by Athanasios Papoulis, example 7-20.
 * @param sig a fpt, the standar deviation of the process.
 * @param d_input a fpt, the input pointer with the samples of the uniform distribution.
 * @param d_output a fpt, the output pointer with the Rayleigh distributed samples.
 * @param n an integer, the number of elements in the input and output pointers.
 */
__global__ void rayleighDistributionKernel(fpt sig,fpt *d_input,fpt *d_output,int n){

    int i = threadIdx.x + blockIdx.x*blockDim.x;

    if(i<n){
      //notice that sqrt and log are from the CUDA math API
      d_output[i] = sig*sqrt(-2*log(d_input[i]));
    }
  }

  template <typename fpt>
/**
 * @brief The kernel riceDistributionkernel() computes samples of the Rice distribution accordingly to the model in "Mobile Radio Channels", by Matthias Patzold, chapter 2.1.2.
 * @param rho a fpt, the amplitude of the Line-of-Sight channel component.
 * @param d_input1 a fpt, the first input pointer with the samples of the Normal distribution.
 * @param d_input2 a fpt, the second input pointer with the samples of the Normal distribution.
 * @param d_output a fpt, the outuput pointer with the Rice distributed samples.
 * @param n an integer, the number of samples in the input and output pointers.
 */
__global__ void riceDistributionKernel(fpt rho, fpt *d_input1, fpt *d_input2, fpt *d_output, int n){

    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i<n) {
      d_output[i] = sqrt(pow(d_input1[i] + rho, 2) + pow(d_input2[i], 2));      
    }
  }

template <typename fpt>
  /**
 * @brief The kernel suzukidistributionkernel() computes the samples of the Suzuki distribution following the model in "Mobile Radio Channels", by Matthias Patzold, chapter 6.1.3, equation (6.54).
 * @param rho a fpt, the amplitude of the Line-of-Sight channel component.
 * @param d_rice a fpt, the input pointer with the samples of the Rice distribution.
 * @param d_log a fpt, the input pointer with the samples of the Lognormal distribution.
 * @param d_output a fpt, the output pointer with the Suzuki distributed samples.
 * @param n an integer, the number of samples in the input and output pointers.
 */  
__global__ void suzukiDistributionKernel(fpt rho, fpt *d_rice, fpt *d_log, fpt *d_output, int n){

    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < n) {
      d_output[i] = d_rice[i] * d_log[i];      
    }    
  }  
#endif
