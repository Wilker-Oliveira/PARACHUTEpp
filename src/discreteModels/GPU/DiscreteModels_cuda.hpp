#ifndef LIBDISCRETEMODELS_CUDA_H
#define LIBDISCRETEMODELS_CUDA_H

//basic libs
#include <cmath>
#include <cstddef>
#include <vector>
#include <algorithm>
//random utilities
#include <random>
#include <chrono>
//Error handling
#include <stdexcept>
//CUDA libs
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
//sublibrary for holding CUDA kernels
#include "DM_kernels.hpp"

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)


/**
 * @brief The namespace gpu_model distincts the GPU implementation from the CPU implementation.
 */

namespace gpu_model{

  template <typename fpt>

  /**
   * @brief The Discreteffchannels class implements the probability density
   * functions for the distributions described in the chapters 2.1 and 6.1 of
   * "Mobile Radio Channels" by Matthias Patzold.
   */
  class DiscreteFFChannels {

  private:
    curandGenerator_t gen; /**< @brief Define the pseudorandom numbers generator.*/
    int ThperBl;

  public:
    DiscreteFFChannels(); /**< @brief The constructor of the class Discreteffchannels.*/
    ~DiscreteFFChannels(); /**< @brief The destructor of the class.*/

    /**
     * @brief Funtion that generates samples of the Gaussian distribution.
     * @param NofSamples a integer, the number of samples to be generated.
     * @param mu a fpt, the mean of the distribution.
     * @param sig a fpt, the standard deviation of the distribution.
     * @return a vector of length Nofsamples with the samples of the Gaussian distribution.  
     */
    std::vector<fpt> gaussDistribution(int NofSamples, fpt mu, fpt sig);

    /**
     * @brief Funtion that generates samples of the Lognormal distribution.
     * @param NofSamples a integer, the number of samples to be generated.
     * @param mu a fpt, the mean of the distribution.
     * @param sig a fpt, the standard deviation of the distribution.
     * @return a vector of length Nofsamples with the samples of the Lognormal distribution.
     */
    std::vector<fpt> lognormDistribution(int NofSamples, fpt mu, fpt sig);

    /**
     * @brief Funtion that generates samples of the Rayleigh distribution.
     * @param NofSamples a integer, the number of samples to be generated.
     * @param sig a fpt, the standard deviation of the distribution.
     * @return a vector of length Nofsamples with the samples of the Rayleigh distribution.
     */
    std::vector<fpt> rayleighDistribution(int NofSamples, fpt sig);

    /**
     * @brief Funtion that generates samples of the Rice distribution.
     * @param NofSamples a integer, the number of samples to be generated.
     * @param sig a fpt, the standard deviation of the distribution.
     * @param rho a fpt, the amplitude of the line-of-sight channel component.
     * @return a vector of length Nofsamples with the samples of the Rice distribution.
     */
    std::vector<fpt> riceDistribution(int NofSamples, fpt sig, fpt rho);
    
    /**
     * @brief Funtion that generates samples of the Suzuki distribution.
     * @param NofSamples a integer, the number of samples to be generated.
     * @param sig0 a fpt, the standard deviation of the Rice distribution.
     * @param rho a fpt, the amplitude of the line-of-sight channel component.
     * @param mu a fpt, the mean of the Lognormal distribution.
     * @param sigu a fpt, the standard deviation of the Lognormal distribution.
     * @return a vector of length Nofsamples with the samples of the Suzuki
     * distribution.
     */  
    std::vector<fpt> suzukiDistribution(int NofSamples, fpt sig0, fpt rho, fpt mu, fpt sigu);
};

  template <typename fpt>
  /**
   * @details
   * Performs the restriction for the fpt template parameter, must be only a Floating Point type.
   * @details
   * Sets the seed and the generator from CURAND for the pseudorandom numbers.
   */
  DiscreteFFChannels<fpt>::DiscreteFFChannels(){

    if(!std::is_floating_point<fpt>()) throw std::runtime_error("fpt should be a Floating Point Type");
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    // Create a pseudorandom number generator. CURAND_RNG_PSEUDO_MTGP32 is a member of the Mersenne Twister family of pseudorandom number generators and has parameters customized for operation on the GPU.
    curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_MTGP32);
    // Set seed.
    curandSetPseudoRandomGeneratorSeed(gen,seed);
    // Set Default number of threads perblock.
    ThperBl = 512;
  }

  template <typename fpt>
  /**
   * @details
   * Performs the destruction of the CURAND pseudorandom numbers generator.
   */  
  DiscreteFFChannels<fpt>::~DiscreteFFChannels(){
      
      curandDestroyGenerator(gen); 
    }

    template <typename fpt>
    /**
     * @details
     * The gaussdistribution() function generates random samples that follow a Gaussian distribution, described by the mean and the standard deviation.
     * @details
     * The PDF of the Gaussian distribution is described by equation (2.30) in "Mobile Radio Channels" by Matthias Patzold.
     */    
   std::vector<fpt> DiscreteFFChannels<fpt>::gaussDistribution(int NofSamples, fpt mu, fpt sig){
    // Creating and zero initializing the return vector.
    std::vector<fpt> ret(NofSamples,0);

    // Creating the data vector pointers.
    fpt *d_gauss, *h_gauss;

    // Allocating the host data vector.
    h_gauss = new fpt[NofSamples]();

    // Allocating the device data vector.
    cudaMalloc(&d_gauss, NofSamples * sizeof(fpt));
    cudaCheckErrors("cudaMalloc failure");    

    // Running a random distribution to generate normally distributed samples.
    if constexpr (std::is_same<fpt, float>()) {curandGenerateNormal(gen, d_gauss, NofSamples, mu, sig);}
    else {curandGenerateNormalDouble(gen, d_gauss, NofSamples, mu, sig);}
    cudaCheckErrors("cudaGenerate failure");
    
    // Copying output from device to host.
    cudaMemcpy(h_gauss, d_gauss, NofSamples * sizeof(fpt), cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy D2H failure");


    // Passing the values of h_ray to ret.
    for(std::size_t i=0;i<ret.size();++i)
      ret[i]=h_gauss[i];
    
    // Freeing the used memory.
    cudaFree(d_gauss);
    delete [] h_gauss;

    return ret;
  }

  template <typename fpt>
  /**
   * @details
   * The lognormdistribution() function generates random samples that follow a Lognormal distribution, described by the mean and the standard deviation in a logarithmic base.
   * @details
   * The PDF of the Lognormal distribution is described by equation (2.51) in "Mobile
Radio Channels" by Matthias Patzold.
   */  
  std::vector<fpt> DiscreteFFChannels<fpt>::lognormDistribution(int NofSamples, fpt mu, fpt sig){
    // Creating and zero initializing the return vector.
    std::vector<fpt> ret(NofSamples, 0);

    // Creating the data vector pointers.
    fpt *d_log, *h_log;

    // Allocating the host data vector.
    h_log = new fpt[NofSamples]();

    // Allocating the device data vector.
    cudaMalloc(&d_log, NofSamples*sizeof(fpt));

    // Running a random distribution to generate uniform distributed samples.    
    if constexpr (std::is_same<fpt, float>()){curandGenerateLogNormal(gen, d_log, NofSamples, mu, sig);}
    else {curandGenerateLogNormalDouble(gen, d_log, NofSamples, mu, sig);}

    // Copying output from device to host.
    cudaMemcpy(h_log, d_log, NofSamples * sizeof(fpt), cudaMemcpyDeviceToHost);

    // Passing the values of h_log to ret.
    for(std::size_t i=0;i<ret.size();++i)
      ret[i]=h_log[i];
    
    // Freeing the used memory.
    cudaFree(d_log);
    delete [] h_log;

    return ret;
        
  }

  template <typename fpt>
  /**
   * @details
   * The rayleighdistribution() function generates random samples that follow a Rayleigh distribution, described by the standard deviation.
   * @details
   * The PDF of the Rayleigh distribution is described by equation (2.40) in "Mobile
Radio Channels" by Matthias Patzold. This implementation uses the Rayleigh distributed random number sequence function in example 7-20 of "Probability, Random Variables and Stochastic Processes" by Athanasios Papoulis, in which the Rayleigh distributed sequence is generated using the Percentile Transformation Method, with the natural logarithm of a uniform random variable.
   */
  std::vector<fpt> DiscreteFFChannels<fpt>::rayleighDistribution(int NofSamples, fpt sig){
    // Creating and zero initializing the return vector.
    std::vector<fpt> ret(NofSamples,0);

    // Creating the data vector pointers.
    fpt *d_uni, *d_ray, *h_ray;

    // Allocating the host data vector.
    h_ray = new fpt[NofSamples]();

    // Allocating the device data vectors.
    cudaMalloc(&d_uni, NofSamples*sizeof(fpt));
    cudaMalloc(&d_ray, NofSamples*sizeof(fpt));  

    // Running a random distribution to generate uniform distributed samples.
    if constexpr(std::is_same<fpt,float>()){ curandGenerateUniform(gen, d_uni, NofSamples);}
    else {curandGenerateUniformDouble(gen, d_uni, NofSamples);}

    // Generating the Rayleigh distributed variables using the uniform samples.
    rayleighDistributionKernel<fpt><<<(NofSamples+(ThperBl-1))/ThperBl,ThperBl>>>(sig, d_uni, d_ray, NofSamples);

    // Copying output from device to host.
    cudaMemcpy(h_ray, d_ray, NofSamples*sizeof(fpt), cudaMemcpyDeviceToHost);

    // Passing the values of h_ray to ret.
    for(std::size_t i=0;i<ret.size();++i)
      ret[i]=h_ray[i];
    
    // Freeing the used memory.
    cudaFree(d_uni);
    cudaFree(d_ray);
    delete [] h_ray;

    return ret;
  }

  template <typename fpt>
  /**
   * @details
   * The ricedistribution() function generates random samples that follow a Rice distribution, described by the standard deviation and the amplitude of the line-of-sight channel component.
   * @details
   * The PDF of the Rice distribution is described by equation (2.44) in "Mobile Radio Channels" by Matthias Patzold.
   */
  std::vector<fpt> DiscreteFFChannels<fpt>::riceDistribution(int NofSamples, fpt sig, fpt rho){
    // Creating and zero initializing the return vector.
    std::vector<fpt> ret(NofSamples, 0);

    // Creating the data vector pointers.
    fpt *d_norm1, *d_norm2, *d_rice, *h_rice;

    // Allocating the host data vector.
    h_rice = new fpt[NofSamples]();

    // Allocating the device data vectors.
    cudaMalloc(&d_norm1, NofSamples*sizeof(fpt));
    cudaMalloc(&d_norm2, NofSamples*sizeof(fpt));    
    cudaMalloc(&d_rice, NofSamples*sizeof(fpt));

    // Running a random distribution to generate uniform distributed samples.
    if constexpr (std::is_same<fpt, float>()){
      curandGenerateNormal(gen, d_norm1, NofSamples, 0.0, sig);
      curandGenerateNormal(gen, d_norm2, NofSamples, 0.0, sig);
    }
    else{
      curandGenerateNormalDouble(gen, d_norm1, NofSamples, 0.0, sig);
      curandGenerateNormalDouble(gen, d_norm2, NofSamples, 0.0, sig);
    }

    // Generating the Rice distributed variables using the normal samples.
    riceDistributionKernel<fpt><<<(NofSamples + (ThperBl - 1)) / ThperBl, ThperBl>>>(rho, d_norm1, d_norm2, d_rice, NofSamples);

    // Copying output from device to host.
    cudaMemcpy(h_rice, d_rice, NofSamples*sizeof(fpt),cudaMemcpyDeviceToHost);

    // Passing the values of h_rice to ret.
    for(std::size_t i=0;i<ret.size();++i)
      ret[i]=h_rice[i];
    
    // Freeing the used memory.
    cudaFree(d_norm1);
    cudaFree(d_norm2);    
    cudaFree(d_rice);
    delete [] h_rice;

    return ret;    
  }

  template <typename fpt>
  /**
   * @details
   * The suzukidistribution() function generates random samples that follow a Suzuki distribution, described by the product of a Rice distributed variable and a Lognormal distributed variable, as shown in equation (6.54) of "Mobile Radio Channels" by Matthias Patzold. The Suzuki distribution is described by the standard deviation of the Rice distribution, the amplitude of the line-of-sight channel component, and the mean and standard deviation of the lognormal distribution.
   * @details
   * The PDF of the Suzuki distribution is described by equation (6.56) in "Mobile Radio Channels" by Matthias Patzold.
   */
  std::vector<fpt>DiscreteFFChannels<fpt>::suzukiDistribution(int NofSamples, fpt sig0, fpt rho, fpt mu, fpt sigu){
    // Creating and zero initializing the return vector.
    std::vector<fpt> ret(NofSamples,0);

    // Creating the data vector pointers.
    fpt *d_rice, *d_log, *d_su, *h_su;

    // Allocating the host data vectors.
    h_su = new fpt[NofSamples]();

    // Allocating the device data vectors.
    cudaMalloc(&d_rice, NofSamples*sizeof(fpt));
    cudaMalloc(&d_log, NofSamples*sizeof(fpt));
    cudaMalloc(&d_su, NofSamples*sizeof(fpt));

    // Generating Rice and Lognormal distributed samples.
    auto rice_vec = DiscreteFFChannels<fpt>::riceDistribution(NofSamples, sig0, rho);
    auto logn_vec = DiscreteFFChannels<fpt>::lognormDistribution(NofSamples, mu, sigu);

    // Copying Rice and Lognormal samples from host to device.
    cudaMemcpy(d_rice, rice_vec.data(), NofSamples*sizeof(fpt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_log, logn_vec.data(), NofSamples*sizeof(fpt), cudaMemcpyHostToDevice);

    // Generating the Suzuki distributed variables using the Rice and Lognormal samples.
    suzukiDistributionKernel<fpt> <<<(NofSamples + (ThperBl - 1)) / ThperBl, ThperBl>>>(rho, d_rice, d_log, d_su, NofSamples);

    // Copying output from device to host.
    cudaMemcpy(h_su, d_su, NofSamples*sizeof(fpt), cudaMemcpyDeviceToHost);

    // Passing the values of h_su to ret.
    for(std::size_t i=0;i<ret.size();++i)
      ret[i]=h_su[i];
    
    // Freeing the used memory.
    cudaFree(d_rice);
    cudaFree(d_log);
    cudaFree(d_su);
    delete [] h_su;
    rice_vec.clear();
    logn_vec.clear();

    return ret;
  }  
}

#endif
