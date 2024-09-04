#ifndef JAKESMODEL_HPP
#define JAKESMODEL_HPP

#ifdef GPU
#include "../Models/GPU/SoSModel_cuda.hpp"
#else
#include "../Models/CPU/SoSModel.hpp"
#endif

#include <iostream>

template<short N, typename fpt>

/**
 * @brief Implementaton of the Jakes Method.
 * @brief This class executes the Jakes method for computing the Doppler Frequencies and Path Gains.
 * @brief This model is described in Section 5.1.5 from Mobile Radio Channels by Matthias Patzold.
 */

class JakesModel: public SoSModel<N,fpt>{

public:

  /**
   * @brief Default constructor for the Jakes power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fmax a float, the maximum Doppler frequency of the process
   * @param i_2 a boolean, true if the process is the complex component of the channel and false otherwise. Default = false
   * @brief Here, the genphases() function is overridden in order to propperly implement equation (5.62) in Section 5.1.5 from Mobile Radio Channels by Matthias Patzold.
   */
  JakesModel(float sig, float fmax, bool i_2 = false){
    if(i_2==false)
      DefineModel(sig, fmax);
    else
      DefineModelImag(sig, fmax);
    genPhases();
  }

  /**
   * @brief Default constructor for the Gauss power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fc a float, the 3dB cut-off frequency
   * @param kc a float, constant to attend the mean power condition
   */
  JakesModel(float sig, float fc, float kc){
    DefineModel(sig,fc,kc);
  }

  /**
   * @brief Defining the function genphases().
   * @brief This function overrides the original genphases() implemented in the SosModel class. Unlike the other models, the Jakes method considers all phases equal to 0 (equation 5.62).
   * @brief Type: void.
   * There are no parameters.
   */
  void genPhases(){
    for(short i=0;i<N;i++) this->phases[i] = 0;
  }

  /**
   * @brief Defining the function DefineModel() for the Jakes PSD of the real component process.
   * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Jakes power spectral density.
   * @brief Type: void.
   * @param sig a float, the standard deviation of the channel
   * @param fmax a float, the maximum Doppler frequency of the channel
   */
  void DefineModel(float sig, float fmax){
 
    float scale=sig/(std::sqrt(N - 0.5));

    for(float n=0;n<N-1;n++)
      this->dopplerFrequencies[n] = fmax*std::cos(M_PI*(n+1.0)/((2*N)-1));
    this->dopplerFrequencies[N-1] = fmax;

    for(float n=0;n<N-1;n++)
      this->pathGains[n] = (2*scale)*std::sin(M_PI*(n+1.0)/(N-1));
    this->pathGains[N-1] = scale;
  }

  /**
   * @brief Defining the function DefineModel() for the Jakes PSD of the imaginary component process.
   * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Jakes power spectral density.
   * @brief Type: void.
   * @param sig a float, the standard deviation of the channel
   * @param fmax a float, the maximum Doppler frequency of the channel
   */
  void DefineModelImag(float sig, float fmax){
    
    float scale=sig/(std::sqrt(N - 0.5));

    for(float n=0;n<N-1;n++)
      this->dopplerFrequencies[n] = fmax*std::cos(M_PI*(n+1.0)/((2*N)-1));
    this->dopplerFrequencies[N-1] = fmax;

    for(float n=0;n<N-1;n++)
      this->pathGains[n] =  (2*scale)* std::cos((M_PI*(n+1.0))/(N-1));
    this->pathGains[N-1] = scale;
  }
  

  
  /**
    * @brief Defining the function DefineModel() for the Gaussian PSD.
    * @brief This function returns a message of error, given that the Jakes method does not perform well for gaussian distribution.
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
  void DefineModel(float sig, float fc, float kc){
    throw std::runtime_error("The process do not follow optimally the Gaussian distribution.");
  }
  
};


#endif
