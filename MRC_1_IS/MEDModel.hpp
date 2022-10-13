#include "SoSModel.hpp"

template<short N, typename fpt>
class MEDModel: public SoSModel<N, fpt> {
  
public:

  void DefineModel(float sig /**< std_dev in lin. */, float fmax){
  for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fmax*(2*n))/(2*N);
  
  for(short n=0;n<N;n++){ 
    this->pathGains[n]= ((2*sig)/sqrt(M_PI))* sqrt(asin(((float)n+1)/N)-asin(((float)n)/N));
  }

}

  void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){
    
    for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fc*kc*(2*n-1))/(2*N);

    for(short n=0;n<N;n++) this->pathGains[n]= (sig*M_SQRT2)*sqrt(erf((n*kc*sqrt(M_LN2))/N) - erf(((n-1)*kc*sqrt(M_LN2))/N));
  }

};

