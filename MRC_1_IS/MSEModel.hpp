#include "SoSModel.hpp"

template<short N, typename fpt>
class MSEModel: public SoSModel<N,fpt>{

public:

  double besselJ0(double x, float step=0.01){
  double I=0,s1=0,s2=0,s3=0;

  //using type deduction to generalize the function f for any datatype.
  auto f=[](auto x,auto o){ return cos(x*cos(o)); };

  for(float i=step;i<(M_PI_2-step);i+=step){
    if(int(i/step)%2!=0) s2+=f(x,i);
    else s3+=f(x,i);
  }
  
  s1=f(x,0)+f(x,M_PI_2);

  I=(s1+4*s2+2*s3)*(step/3);
  
  return M_2_PI*I;
};


  double JakesIntegral(fpt dfreq, float tmax, float fmax, float step =0.001){


    double I=0,s1=0,s2=0,s3=0;
    auto f=[this,fmax](fpt dfreq, auto tau){return besselJ0(2*M_PI*fmax*tau)*cos(2*M_PI*dfreq*tau);};
  
    for(float i=step;i<(tmax-step);i+=step){
      if(int(i/step)%2!=0) s2+=f(dfreq,i);
      else s3+=f(dfreq,i);
    }
  
    s1=f(dfreq,0)+f(dfreq,tmax);
    
    I=(s1+4*s2+2*s3)*(step/3);
    
    return I;
}

  double GaussIntegral(fpt dfreq, float tmax, float fc, float step =0.001){

    double I=0,s1=0,s2=0,s3=0;
    auto e=[fc](auto tau){return exp(-pow((M_PI*fc*tau),2)/(M_LN2));};

    auto f=[e](fpt dfreq, auto tau){return e(tau)*cos(2*M_PI*dfreq*tau);};

    for(float i=step;i<(tmax-step);i+=step){
      if(int(i/step)%2!=0) s2+=f(dfreq,i);
      else s3+=f(dfreq,i);
    }

    s1=f(dfreq,0)+f(dfreq,tmax);

    I=(s1+4*s2+2*s3)*(step/3);

    return I;
  }

    void DefineModel(float sig /**< std_dev in lin. */, float fmax){
    float tmax = N/(2*fmax);

    for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fmax*(2*(n+1)-1))/(2*N);

    for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax) * JakesIntegral(this->dopplerFrequencies[n],tmax,fmax));
    }

    void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){

      float tmax = N/(2*kc*fc);

      for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fc*kc*(2*(n+1)-1))/(2*N);

      for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax)* GaussIntegral(this->dopplerFrequencies[n],tmax,fc));
    }

};
