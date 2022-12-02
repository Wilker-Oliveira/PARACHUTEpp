#ifndef MSEMODEL_HPP
#define MSEMODEL_HPP

#include "SoSModel.hpp"

template<short N, typename fpt>

/**
 * @brief Implementaton of the Mean Square Error Method.
 * @brief This class executes the MSE method for computing the Doppler Frequencies and Path Gains.
 * @brief This model is described in Section 5.1.2 from Mobile Radio Channels by Matthias Patzold.
 */

class MSEModel: public SoSModel<N,fpt>{

public:
  
  /**
   * @brief Default constructor for the Jakes power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fmax a float, the maximum Doppler frequency of the process
   */
  MSEModel(float sig, float fmax){
    DefineModel(sig, fmax);
    this->genPhases();
  }
  
  /**
   * @brief Default constructor for the Gauss power spectral density form.
   * @param sig a float, standard deviation of the process
   * @param fc a float, the 3dB cut-off frequency
   * @param kc a float, constant to attend the mean power condition
   */
  MSEModel(float sig, float fc, float kc){
    DefineModel(sig, fc, kc);
    this->genPhases();
  }
  
  /**
   * @brief Defining the function besselJ0().
   * @brief This function implements the Bessel function of the zero-th order J0(x) following the equation (4.4) from Introduction to Bessel Functions by Frank Bowman, that expresses J0(x) as a definite integral.
   * @brief The integration is performed numerically through the Simpson's 1/3 rule.
   * @brief Type: double.
   * @param x a double, argument of the Bessel function of the zero-th order
   * @param step a float, step of the numerical integration, default value 0.01
   * @return The value of the Bessel function J0(x).
   */
  double besselJ0(double x, float step=0.01) const {
  
    double I=0,s1=0,s2=0,s3=0;
    /** Using type deduction to generalize the function f for any datatype. */
    auto f=[](auto x,auto o){ return cos(x*cos(o)); };

    for(float i=step;i<(M_PI_2-step);i+=step){
      if(int(i/step)%2!=0) s2+=f(x,i);
      else s3+=f(x,i);
    }
  
    s1=f(x,0)+f(x,M_PI_2);

    I=(s1+4*s2+2*s3)*(step/3);
    
    return M_2_PI*I;
  };

  /**
   * @brief Defining the function JakesIntegral()
   *
   * @brief This function implements the integral in equation (5.22) from Mobile Radio Channels.
   * @brief The integration is performed numerically through the Simpson's 1/3 rule.
   * @brief Type: double.
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param tmax a float, represents the appropriate time-lag interval and upper limit of the integral
   * @param fmax a float, the maximum Doppler frequency of the process
   * @param step a float, step of the numerical integration, default value 0.001
   * @return the integral value
   */
  double JakesIntegral(fpt dfreq, float tmax, float fmax, float step =0.001) const {


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

  /**
   * @brief Defining the function GaussianIntegral()
   *
   * @brief Implements the integral in equation (5.24) from Mobile Radio Channels.
   * @brief The integration is performed numerically through the Simpson's 1/3 rule.
   * @brief Type: double.
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param tmax a float, the appropriate time-lag interval and upper limit of the integral
   * @param fc a float, the 3dB cut-off frequency
   * @param step a float, step of the numerical integration, default value 0.001
   * @return the integral value
   */
  double GaussIntegral(fpt dfreq, float tmax, float fc, float step =0.001) const {

    double I=0,s1=0,s2=0,s3=0;
    auto e=[fc](auto tau){return exp(-pow((M_PI*fc*tau),2)/(M_LN2));};

    auto f=[&e](fpt dfreq, auto tau){return e(tau)*cos(2*M_PI*dfreq*tau);};

    for(float i=step;i<(tmax-step);i+=step){
      if(int(i/step)%2!=0) s2+=f(dfreq,i);
      else s3+=f(dfreq,i);
    }

    s1=f(dfreq,0)+f(dfreq,tmax);

    I=(s1+4*s2+2*s3)*(step/3);

    return I;
  }

   /**
    * @brief Defining the function DefineModel() for the Jakes PSD.
    * @brief This function computes the doppler frequencies and path gains of the MSE method applied on the Jakes power spectral density.
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel in linear scale
    * @param fmax a float, the maximum Doppler frequency of the channel
    */
    void DefineModel(float sig, float fmax){
    float tmax = N/(2*fmax);

    for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fmax*(2*(n+1)-1))/(2*N);

    for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax) * JakesIntegral(this->dopplerFrequencies[n],tmax,fmax));
    }

   /**
    * @brief Defining the function DefineModel() for the Gaussian PSD.
    * @brief This function computes the doppler frequencies and path gains of the MSE method applied on the Gaussian power spectral density
    * @brief Type: void.
    * @param sig a float, the standard deviation of the channel in linear scale
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
    void DefineModel(float sig, float fc, float kc){

      float tmax = N/(2*kc*fc);

      for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fc*kc*(2*(n+1)-1))/(2*N);

      for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax)* GaussIntegral(this->dopplerFrequencies[n],tmax,fc));
    }

};

#endif
