#ifndef MEDSMODEL_HPP
#define MEDSMODEL_HPP

#include "SoSModel.hpp"

template<short N, typename fpt>

class MEDSModel: public SoSModel<N,fpt>{

protected:

  fpt GaussianPSDe(fpt dfreq, fpt n, float fc){
    return ( ((2*(n+1)-1)/((fpt)2*N))-std::erf(dfreq*(std::sqrt( std::numbers::ln2_v<fpt> )/fc) ) );
  }

public:

  MEDSModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }
  
  MEDSModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }

  template<typename F>
  
  /**
   * @brief Declaring the function bisecmethod(). 
   */
  fpt bisecMethod(fpt a, fpt b, fpt tol/*tolerance*/, F f, short lim=100);


  void DefineModel(float sig, float fmax){
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));
    for(short n=0;n<N;n++) this->dopplerFrequencies[n]=fmax*std::sin((M_PI/(2*N))*(n-0.5));
  }
  
  void DefineModel(float sig, float fc, float kc){
    //Beta from MRC 3.68
    fpt aux=0, beta=2*pow(M_PI*fc*sig,2)/log(2);
    
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));

    for(short n=0;n<N-1;n++){      
      auto GaussF= [n,fc,this](fpt dfreq){return this->GaussianPSDe(dfreq, n, fc);};
      this->dopplerFrequencies[n]= this->bisecMethod(0, 3*fc/std::sqrt(std::numbers::ln2_v<float>), 0.005, GaussF);
      aux+=pow(this->dopplerFrequencies[n],2);
    }

    this->dopplerFrequencies[N-1] = sqrt((beta*N)/pow(2*M_PI*sig,2) - aux);

  }


};

template<short N, typename fpt>
template<typename F>

/**
 * @brief Defining the function bisecmethod().
 * @brief This function implements the bisection method to find the roots of a given expression in a determined interval.
 * @brief The parameters of this function must attend the bisection method criteria.
 * @param a a fpt, the inferior limit of the interval
 * @param b a fpt, the superior limit of the interval
 * @param tol a fpt, the stopping criterion of the method
 * @param f as F, a typename template, the function for which the root must be found
 * @param lim a short, the maximum number of iteractions to be executed in order to guarantee convergence
 * @return The root of the function in the interval [a,b].
 */
fpt MEDSModel<N, fpt>::bisecMethod(fpt a, fpt b, fpt tol, F f, short lim){
  
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

  return res;
};



#endif
