#include "SoSModel.hpp"
#include "../1_Rayleigh_and_Rice/SpectralProperties.hpp"

template<short N, typename fpt>
class MSEModel: public SoSModel<N,fpt>{

public:

  double JakesIntegral(fpt dfreq, float tmax, float fmax, float step =0.01){


    double I=0,s1=0,s2=0,s3=0;
    auto f=[fmax](fpt dfreq, auto tau){return besselJ0(2*M_PI*fmax*tau)*cos(2*M_PI*dfreq*tau);};
  
    for(float i=step;i<(tmax-step);i+=step){
      if(int(i/step)%2!=0) s2+=fi(dfreq,i);
      else s3+=f(dfreq,i);
    }
  
    s1=f(dfreq,0)+ft(dfreq,tmax);
    
    I=(s1+4*s2+2*s3)*(step/3);
  
    return I;
}

  double GaussIntegral(fpt dfreq, float tmax, float fc, float kc, float step =0.01);

    void DefineModel(float sig /**< std_dev in lin. */, float fmax){
    float tmax = N/(2*fmax);

    for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fmax*(2*n))/(2*N);

     for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax)* JakesIntegral(this->dopplerFrequencies[n],tmax,fmax));

    }

};
