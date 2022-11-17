#include "MEAModel.hpp"

template<short N, typename fpt>

class MCModel: public MEAModel<N, fpt>{

private:
  
  /* Uniform random variable*/
  fpt u_rv(){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<fpt> distribution(0.0, 1.0);

    return distribution(generator);
  }
  
protected:
  
  fpt gaussianPSDe(fpt dfreq, float fc){
    return u_rv() - std::erf((dfreq/fc)*std::sqrt( std::numbers::ln2_v<fpt> ));
  }
  
public:

  MCModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }

  MCModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }

  void DefineModel(float sig, float fmax){
    for(short n=0;n<N;n++) this->pathGains[n] = sig*std::sqrt(2/((float)N)); 
    for(short n=0;n<N;n++) this->dopplerFrequencies[n] = fmax*std::sin(M_PI_2*u_rv());
  }

  void DefineModel(float sig, float fc, float kc){
    for(short n=0;n<N;n++) this->pathGains[n] = sig*std::sqrt(2/((float)N)); 
    for(short n=0;n<N;n++){
      auto GaussF = [fc,this](fpt dfreq){return this->gaussianPSDe(dfreq, fc);};

      this->dopplerFrequencies[n] = this->bissecMethod(0, 3*fc/std::sqrt(std::numbers::ln2_v<float>), 0.005, GaussF);

    }
  }
};
