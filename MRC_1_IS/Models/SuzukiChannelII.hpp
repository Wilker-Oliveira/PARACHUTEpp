#ifndef SUZUKICHANNELII_HPP
#define SUZUKICHANNELII_HPP

#include "MEDSModel.hpp"

template<short N1, short N3, typename fpt>

/**
 * @brief Class implementing the extended Suzuki process of type II.
 * This model is described in Section 6.2 from Mobile Radio Channels by Matthias Patzold. 
 * The MEDS method is employed for computing the path gains and Doppler Frequencies.
 */

class SuzukiChannelII{

protected:
  /**Declaring the Rayleigh process. */
  MEDSModel<N1, fpt> *u1;
  
  /** Declaring the Gaussian process. */
  MEDSModel<N3, fpt> *u3;
  
  fpt p; /**< The amplitude of the LoS component.*/
  fpt fp; /**< The Doppler frequency of the LoS component.*/
  fpt theta_p; /**< The phase shift of the LoS component.*/

  fpt sig3; /**< The fitting standard deviation of the Gaussian process.*/
  fpt m3; /**< The fitting mean value of the Gaussian process.*/

  fpt theta_o; /**< The phase shift between the processes \f$\tilde{\mu}_1\f$ and \f$\tilde{\mu}_2\f$. */

public:

  /**
   * @brief Default constructor of the Suzuki channel of type I.
   * @param sig a float, the standard deviation of the Rayleigh processes.
   * @param fmax a float, the maximum Doppler frequency ot the Rayleigh processes.
   * @param ko a float, the ratio of the minimum and maximum doppler frequencies, with values in the interval \f$[0,1]\f$.
   * @param theta_o a fpt, the phase shift between the two Rayleigh processes.
   * @param sig3 a fpt, the fitting standard deviation of the Gaussian process.
   * @param fc a float, the 3dB cut-off frequency of the Gaussian process.
   * @param kc a float, constant to attend the mean power condition.
   * @param m3 a fpt, the fitting mean of the Gaussian process.
   * @param p a fpt, the amplitude of the line-of-sight component.
   * @param fp a fpt, the Doppler frequency of the line-of-sight component.
   * @param theta_p a fpt, the phase shift of the line-of-sight component.
   */
  SuzukiChannelII(float sig, float fmax, float ko, fpt theta_o, fpt sig3, float fc, float kc, fpt m3, fpt p, fpt fp, fpt theta_p){
    
    u1 = new MEDSModel<N1, fpt>(sig,fmax,true,ko);
    u3 = new MEDSModel<N3, fpt>(1,fc, kc);
    this->p=p;
    this->fp=fp;
    this->theta_p=theta_p;
    this->sig3=sig3;
    this->m3=m3;
    this->theta_o=theta_o;
  }

  /**
   * @brief Declaring the CalChannel() function.
   * @brief This function computes the real, imaginary and lognormal processes that compose the Suzuki channel of type I, 
   * accordingly with equation (6.54) from Mobile Radio Channels by Matthias Patzold.
   * @brief Type: vector of type fpt.
   * @param time a vector of floats, the considered time interval of the process.
   * @return a vector containing the Suzuki process of type I values for each instant of time.
   */
  std::vector<fpt> CalcChannel(const std::vector<float> &time) const {
    fpt re_sq=0,im_sq=0, lgn=0;
    std::vector<fpt> res;
    
    std::vector<fpt> u1_cos = u1->CalcProcess(time);
    u1->addPhase(-theta_o);
    std::vector<fpt> u2_sin = u1->CalcProcess(time);

    std::vector<fpt> u3_cos = u3->CalcProcess(time);

    for(int i=0;i<time.size();i++){
      re_sq = std::pow(u1_cos[i] + p*std::cos(theta_p),2);
      im_sq = std::pow(u2_sin[i] + p*std::sin(theta_p),2);

      lgn=std::exp(u3_cos[i]*sig3+m3);

      res.push_back((std::sqrt(re_sq+im_sq))*lgn);
    }

    return res;
  }
  
};

#endif
