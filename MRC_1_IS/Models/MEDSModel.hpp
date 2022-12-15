#ifndef MEDSMODEL_HPP
#define MEDSMODEL_HPP

#include "SoSModel.hpp"

template<short N, typename fpt>

/**
 * @brief Implementaton of the Method of Exact Doppler Spread.
 * @brief This class executes the MEDS method for computing the Doppler Frequencies and Path Gains.
 * @brief This model is described in Section 5.1.7 from Mobile Radio Channels by Matthias Patzold.
 */

class MEDSModel: public SoSModel<N,fpt>{

protected:

  /**
   * @brief Template function for gaussian power spectral density equation.
   * @brief This function is the implementaton of equation (5.89a) in Mobile Radio Channels by Matthias Patzold.
   * It is used to determine the doppler frequencies of the process.
   * @brief Type: fpt.
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param n a fpt, the index of the frequency
   * @param fc a float, the 3dB cut-off frequency
   * @return The value of expression.
   */
  fpt GaussianPSDe(fpt dfreq, fpt n, float fc){
    return ( ((2*(n+1)-1)/((fpt)2*N))-std::erf(dfreq*(std::sqrt( std::numbers::ln2_v<fpt> )/fc) ) );
  }
  fpt ko;

public:

  /**
   * @brief Default constructor for the Jakes power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fmax a float, the maximum Doppler frequency of the process
   * @param halfpower a boolean, true when simulating a Suzuki channel, false otherwise; default = false
   * @param ko a float, the ratio of the minimum and maximum doppler frequencies, with values in the interval \f$[0,1]\f$; default = 1
   */
  MEDSModel(float sig, float fmax, bool halfpower=false,float ko=1){
    this->ko=ko;
    DefineModel(sig, fmax);
    if(halfpower==true) ScalePathGains(1/M_SQRT2);
    this->genPhases();
  }

  /**
   * @brief Default constructor for the Gauss power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fc a float, the 3dB cut-off frequency
   * @param kc a float, constant to attend the mean power condition
   */
  MEDSModel(float sig, float fc, float kc){
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
   * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Jakes power spectral density.
   * @brief Type: void.
   * @param sig a float, the standard deviation of the channel
   * @param fmax a float, the maximum Doppler frequency of the channel
   */
  void DefineModel(float sig, float fmax){
    float Nscale=N/(M_2_PI*asin(this->ko));
    for(short n=0;n<N;n++) this->pathGains[n]=sig*std::sqrt(2/((float)Nscale));
    for(short n=0;n<N;n++) this->dopplerFrequencies[n]=fmax*std::sin((M_PI/(2*Nscale))*(n-0.5));
  }

  /**
    * @brief Defining the function DefineModel() for the Gaussian PSD.
    * @brief This function computes the doppler frequencies and path gains of the MCM method applied on the Gaussian power spectral density
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
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

  /**
   * @brief Defining the ScalePathGains() function.
   * @brief This function modifies the path gains computed to attend the power adaptation of the Suzuki process.
   * @brief Type: void.
   * @param scaleFactor a float, the component that performs the scaling of the path gains.
   */
  void ScalePathGains(float scaleFactor){
    for(short n=0;n<N;n++) this->pathGains[n]*=scaleFactor;
  }
  /**
   * @brief Defining the addPhase() function.
   * @brief This function adds a phase shift on the process.
   * @brief Type: void.
   * @param theta a float, the considered phase shift.
   */
  void addPhase(float theta){
    for(short n=0;n<N;n++) this->phases[n]+=theta;
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
