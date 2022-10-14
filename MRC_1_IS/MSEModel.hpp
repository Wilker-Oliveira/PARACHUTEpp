#include "SoSModel.hpp"
#include "../1_Rayleigh_and_Rice/SpectralProperties.hpp"

template<short N, typename fpt>

double SoSModel<N,fpt>::integral(typename fpt dfreq, float tmax, float fmax, step =0.01){

    auto f=[](auto tau){return besselJ0(2*M_PI*fmax*tau)*cos(2*M_PI*dfreq*tau);};

    for(float i=step;i<(tmax-step);i+=step){
        if(int(i/step)%2!=0) s2+=f(i);
        else s3+=f(i);
    }

    s1=f(0)+f(tmax);

    I=(s1+4*s2+2*s3)*(step/3);

    return I;
}


template<short N, typename fpt>
class MSEModel: public SoSModel<N,fpt>{

public:

    void DefineModel(float sig /**< std_dev in lin. */, float fmax){
    float tmax = N/(2*fmax);

    for(short n=0;n<N;n++) this->dopplerFrequencies[n]= (fmax*(2*n))/(2*N);

     for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax)* integral(dopplerFrequencies[n],tmax,fmax));

    }

};
