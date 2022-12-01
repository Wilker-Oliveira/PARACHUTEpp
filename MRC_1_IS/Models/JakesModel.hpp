#include "SoSModel.hpp"

template<short N, typename fpt>

// Document me

class JakesModel: public SoSModel<N,fpt>{

public:
  
  JakesModel(float sig, float fmax, bool i_2 = false){
    DefineModel(sig, fmax, i_2);
    genPhases();
  }

  JakesModel(float sig, float fc, float kc){
    DefineModel(sig,fc,kc);
  }

  void genPhases(){
    for(short i=0;i<N;i++) this->phases[i] = 0;
  }

  void DefineModel(float sig, float fmax, bool i_2 = false){
    
    for(short n=0;n<N;n++){
      if(n==N-1) this->dopplerFrequencies[n] = fmax*std::cos(M_PI*(n+1)/((2*N)-1));
      else this->dopplerFrequencies[n] = fmax;
    }

    if(i_2==false)
      for(short n=0;n<N-1;n++)
	this->pathGains[n] = (2*sig/(std::sqrt(N - 0.5)))*std::sin(M_PI*(n+1)/(N-1));
      
    else
      for(short n=0;n<N-1;n++)
	this->pathGains[n] =  (2*sig/(std::sqrt(N - 0.5)))*std::cos(M_PI*(n+1)/(N-1));
    
    
    this->pathGains[N-1] = sig/(std::sqrt(N - 0.5));
  }

  void DefineModel(float sig, float fc, float kc){
    throw std::runtime_error("The process do not follow optimally the Gaussian distribution.");
  }
  
};
