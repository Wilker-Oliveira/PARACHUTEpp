#include "SoSModel.hpp"

template<short N, typename fpt>

/**
 * \brief Implementaton of the Method of Equal Distances.
 *
 * This class executes the MED method for computing the Doppler Frequencies and Path Gains.
 * This model is described in Section 5.1.1 from Mobile Radio Channels by Matthias Patzold.
 */

class MEDModel: public SoSModel<N, fpt> {
  
public:

 /**
  * \brief Defining the function DefineModel() for the Jakes PSD
  *
  * Computes the doppler frequencies and path gains of the MSE method applied on the Jakes power spectral density
  * type: void
  * @param sig a float, the standard deviation of the channel
  * @param fmax a float, the maximum Doppler frequency of the channel
  */
  void DefineModel(float sig /**< std_dev in lin. */, float fmax){
  for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fmax*(2*(n+1)-1))/(2*N);
  
  for(short n=0;n<N;n++){ 
    this->pathGains[n]= ((2*sig)/sqrt(M_PI))* sqrt(asin(((float)n+1)/N)-asin(((float)n)/N));
  }

}

 /**
  * \brief Defining the function DefineModel() for the Gaussian PSD
  *
  * Computes the doppler frequencies and path gains of the MSE method applied on the Gaussian power spectral density
  * type: void
  * @param sig a float, the standard deviation of the channel
  * @param fc a float, the 3dB cut-off frequency
  * @param kc a float, constant to attend the mean power condition
  */
  void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){
    
    for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fc*kc*(2*(n+1)-1))/(2*N);

    for(short n=0;n<N;n++) this->pathGains[n]= (sig*M_SQRT2)*sqrt(erf(((n+1)*kc*sqrt(M_LN2))/N) - erf(((n)*kc*sqrt(M_LN2))/N));
  }

};

