#ifndef MEAMODEL_HPP
#define MEAMODEL_HPP

#include "SoSModel.hpp"

template<short N, typename fpt>

/**
 * @brief Implementaton of the Method of Equal Areas.
 * @brief This class executes the MEA method for computing the Doppler Frequencies and Path Gains.
 * @brief This model is described in Section 5.1.3 from Mobile Radio Channels by Matthias Patzold.
 */

class MEAModel: public SoSModel<N,fpt>{

protected: 

  /**
   * @brief Template function for gaussian power spectral density equation.
   * @brief This function is the implementaton of equation (5.42) in Mobile Radio Channels by Matthias Patzold.
   * It is used to determine the doppler frequencies of the process.
   * @brief Type: fpt.
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param n a fpt, the index of the n-th sinusoid of the process
   * @param fc a float, the 3dB cut-off frequency
   * @return The value of expression.
   */
  fpt GaussianPSDe(fpt dfreq, fpt n, float fc){
    return ( ((n+1)/((fpt)N))-std::erf(dfreq*(std::sqrt( std::numbers::ln2_v<fpt> )/fc) ) );
  }

public:

  /**
   * @brief Default constructor for the Jakes power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fmax a float, the maximum Doppler frequency of the process
   */
  MEAModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }

  /**
   * @brief Default constructor for the Gauss power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fc a float, the 3dB cut-off frequency
   * @param kc a float, constant to attend the mean power condition
   */  
  MEAModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }

  template<typename F>
  
  /**
   * @brief Declaring the function bisecmethod(). 
   */
  fpt bisecMethod(fpt a, fpt b, fpt tol/*tolerance*/, F f, short lim=100);

  /**
   * @brief Defining the function DefineModel() for the Jakes PSD.
   * @brief This function computes the doppler frequencies and path gains of the MEA method applied on the Jakes power spectral density.
   * @brief Type: void.
   * @param sig a float, the standard deviation of the channel in linear scale
   * @param fmax a float, the maximum Doppler frequency of the channel
   */
  void DefineModel(float sig, float fmax){
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));

    for(short n=0;n<N;n++) this->dopplerFrequencies[n]=fmax*std::sin((M_PI*(n+1))/(2*((float)N)));
    
}

  /**
    * @brief Defining the function DefineModel() for the Gaussian PSD.
    * @brief This function computes the doppler frequencies and path gains of the MEA method applied on the Gaussian power spectral density
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel in linear scale
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
  void DefineModel(float sig, float fc, float kc){
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)N));

    for(short n=0;n<N;n++){      
      auto GaussF= [n,fc,this](fpt dfreq){return this->GaussianPSDe(dfreq, n, fc);};

      this->dopplerFrequencies[n]= this->bisecMethod(0, 3*fc/std::sqrt(std::numbers::ln2_v<float>), 0.005, GaussF);
    }
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
fpt MEAModel<N, fpt>::bisecMethod(fpt a, fpt b, fpt tol, F f, short lim){
  
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
