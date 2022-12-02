#ifndef MCMODEL_HPP
#define MCMODEL_HPP

#include "MEAModel.hpp"

template<short N, typename fpt>

/**
 * @brief Implementaton of the Monte Carlo Method.
 * @brief This class executes the MCM method for computing the Doppler Frequencies and Path Gains.
 * @brief This model is described in Section 5.1.4 from Mobile Radio Channels by Matthias Patzold.
 */

class MCModel: public MEAModel<N, fpt>{

private:
  
  /**
   * @brief Defining the function u_rv().
   * @brief This function computes a pseudorandom number uniformely distributed in the interval [0.0, 1.0].
   * It is used in the functions that determine the Doppler frequencies of the process in the MCM model.
   * @brief Type: fpt.
   * @brief There are no parameters.
   * @return A sample of the uniform distribution in the interval [0.0, 1.0].
   */
  fpt u_rv(){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<fpt> distribution(0.0, 1.0);

    return distribution(generator);
  }
  
protected:

  /**
   * @brief Template function for gaussian power spectral density equation.
   * @brief This function is the implementaton of equation (5.59) in Mobile Radio Channels by Matthias Patzold.
   * It is used to determine the doppler frequencies of the process.
   * @brief Type: fpt.
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param fc a float, the 3dB cut-off frequency
   * @return The value of expression.
   */
  fpt gaussianPSDe(fpt dfreq, float fc){
    return u_rv() - std::erf((dfreq/fc)*std::sqrt( std::numbers::ln2_v<fpt> ));
  }
  
public:

  /**
   * @brief Default constructor for the Jakes power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fmax a float, the maximum Doppler frequency of the process
   */
  MCModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }

  /**
   * @brief Default constructor for the Gauss power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fc a float, the 3dB cut-off frequency
   * @param kc a float, constant to attend the mean power condition
   */
  MCModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }

  /**
   * @brief Defining the function DefineModel() for the Jakes PSD.
   * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Jakes power spectral density.
   * @brief Type: void.
   * @param sig a float, the standard deviation of the channel
   * @param fmax a float, the maximum Doppler frequency of the channel
   */
  void DefineModel(float sig, float fmax){
    for(short n=0;n<N;n++) this->pathGains[n] = sig*std::sqrt(2/((float)N)); 
    for(short n=0;n<N;n++) this->dopplerFrequencies[n] = fmax*std::sin(M_PI_2*u_rv());
  }

  /**
    * @brief Defining the function DefineModel() for the Gaussian PSD.
    * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Gaussian power spectral density
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
  void DefineModel(float sig, float fc, float kc){
    for(short n=0;n<N;n++) this->pathGains[n] = sig*std::sqrt(2/((float)N)); 
    for(short n=0;n<N;n++){
      auto GaussF = [fc,this](fpt dfreq){return this->gaussianPSDe(dfreq, fc);};

      this->dopplerFrequencies[n] = this->bissecMethod(0, 3*fc/std::sqrt(std::numbers::ln2_v<float>), 0.005, GaussF);

    }
  }
};

#endif
