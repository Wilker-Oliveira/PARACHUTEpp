#include "SoSModel.hpp"

template<short N, typename fpt>

/**
 * \brief Implementaton of the Mean Square Error Method.
 *
 * This class executes the MSE method for computing the Doppler Frequencies and Path Gains.
 * This model is described in Section 5.1.2 from Mobile Radio Channels by Matthias Patzold.
 */

class MSEModel: public SoSModel<N,fpt>{

public:

  /**
   * \brief Defining the function besselJ0()
   *
   * Implements the Bessel function of the zero-th order J0(x) following the equation (4.4) from Introduction to Bessel Functions by Frank Bowman, that expresses J0(x) as a definite integral.
   * Performs a numerical integration through the Simpson's 1/3 rule
   * type: double
   * @param x a double, argument of the Bessel function of the zero-th order
   * @param step a float, step of the numerical integration, default value 0.01
   * @return the value of the Bessel function J0(x)
   */
  double besselJ0(double x, float step=0.01){
  double I=0,s1=0,s2=0,s3=0;

  /** using type deduction to generalize the function f for any datatype. */
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
   * \brief Defining the function JakesIntegral()
   *
   * Implements the integral in equation (5.22) from Mobile Radio Channels.
   * Performs a numerical integration through the Simpson's 1/3 rule
   * type: double
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param step a float, the appropriate time-lag interval and upper limit of the integral
   * @param fmax a float, the maximum Doppler frequency of the channel
   * @param step a float, step of the numerical integration, default value 0.001
   * @return I the integral value
   */
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

  /**
   * \brief Defining the function GaussianIntegral()
   *
   * Implements the integral in equation (5.24) from Mobile Radio Channels.
   * Performs a numerical integration through the Simpson's 1/3 rule
   * type: double
   * @param dfreq a fpt, represents the dopplerFrequencies elements
   * @param step a float, the appropriate time-lag interval and upper limit of the integral
   * @param fc a float, the 3dB cut-off frequency
   * @param step a float, step of the numerical integration, default value 0.001
   * @return I the integral value
   */
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

   /**
    * \brief Defining the function DefineModel() for the Jakes PSD
    *
    * Computes the doppler frequencies and path gains of the MSE method applied on the Jakes power spectral density
    * type: void
    * @param sig a float, the standard deviation of the channel
    * @param fmax a float, the maximum Doppler frequency of the channel
    */
    void DefineModel(float sig /**< std_dev in lin. */, float fmax){
    float tmax = N/(2*fmax);

    for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fmax*(2*(n+1)-1))/(2*N);

    for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax) * JakesIntegral(this->dopplerFrequencies[n],tmax,fmax));
    }

   /**
    * \brief Defining the function DefineModel() for the Gaussian PSD
    *
    * Computes the doppler frequencies and path gains of the MSE method applied on the Gaussian power spectral density
    * type: void
    * @param sig a float, the standard deviation of the channel
    * @param fc a float, the 3dB cut-off frequency
    * @param kc a float, constant to attend the mean power condition
    */
    void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){

      float tmax = N/(2*kc*fc);

      for(short n=0;n<N;n++) this->dopplerFrequencies[n] = (fc*kc*(2*(n+1)-1))/(2*N);

      for(short n=0;n<N;n++) this->pathGains[n] = (2*sig)*sqrt( (1.0/tmax)* GaussIntegral(this->dopplerFrequencies[n],tmax,fc));
    }

};
