#include "SoSModel.hpp"
//#include <iostream>
template<short N, typename fpt>

class MEAModel: public SoSModel<N,fpt>{

protected: 

  //template function for gaussian PSD equation
  fpt GaussianPSDe(fpt dfreq, fpt n, float fc){
    return ( (n/((fpt)N))-std::erf(dfreq*(std::sqrt( std::numbers::ln2_v<fpt> )/fc) ) );
  }

public:


  MEAModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }

  MEAModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }

  //bissec method receiving the cost function
  //change a,b to be a tuble container
  template<typename F>
  fpt bissecMethod(fpt a, fpt b, fpt tol/*tolerance*/, F f, short lim=100);

  //DefineModel for jakes PSD
  void DefineModel(float sig /**< std_dev in lin. */, float fmax){
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));

    for(short n=0;n<N;n++) this->dopplerFrequencies[n]=fmax*std::sin((M_PI*n)/(2*((float)N)));
    
}

  //DefineModel for Gaussian PSD
  void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));

    for(short n=0;n<N;n++){      
      auto GaussF= [n,fc,this](fpt dfreq){return this->GaussianPSDe(dfreq, n, fc);};

      this->dopplerFrequencies[n]= this->bissecMethod(0, 3*fc/std::sqrt(std::numbers::ln2_v<float>), 0.005, GaussF);
    }
  }

};

template<short N, typename fpt>
template<typename F>
fpt MEAModel<N, fpt>::bissecMethod(fpt a, fpt b, fpt tol, F f, short lim){
  
  fpt res=(a+b)/2;
  short c=0;
  
  if(fabs(f(a))<=tol) return a;

  if(fabs(f(b))<=tol) return b;

  while(fabs(f(res))>tol && c<lim){
    if(f(a)*f(res)<0) b=res;
    else a=res;

    res=(a+b)/2;
    c++;
  }

  //std::cout<<"c ->"<<c<<'\n';
  //std::cout<<"x ->"<<res<<'\n';
  //std::cout<<"f(x) ->"<<f(res)<<'\n';
  
  return res;
};
