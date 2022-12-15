#ifndef SUZUKICHANNELI_HPP
#define SUZUKICHANNELI_HPP

#include "MEDSModel.hpp"

template<short N1,short N2, short N3, typename fpt>

class SuzukiChannelI{

protected:
  //Rayleigh process
  MEDSModel<N1, fpt> *u1;
  MEDSModel<N2, fpt> *u2;
  //Gaussian process
  MEDSModel<N3, fpt> *u3;
  
  //amplitude, dopplerfrequency and phase shift of the line-of-sight component
  fpt p, fp, theta_p; 

  fpt sig3,m3;
public:

  SuzukiChannelI(float sig, float fmax, float ko, fpt sig3,float fc, float kc, float m3, fpt p, fpt fp, fpt theta_p){
    u1 = new MEDSModel<N1, fpt>(sig,fmax,true);
    u2 = new MEDSModel<N2, fpt>(sig,fmax,true,ko);
    u3 = new MEDSModel<N3, fpt>(1,fc, kc);
    this->p=p;
    this->fp=fp;
    this->theta_p=theta_p;
    this->sig3=sig3;
    this->m3=m3;
  }

  std::vector<fpt> CalcChannel(const std::vector<float> &time) const {
    fpt re_sq=0,im_sq=0, lgn=0;
    std::vector<fpt> res;
    
    std::vector<fpt> u1_cos = u1->CalcProcess(time);
    std::vector<fpt> u2_cos = u2->CalcProcess(time);
    std::vector<fpt> u1_sin = u1->CalcHilbertProcess(time);
    u2->ScalePathGains(-1.0);
    std::vector<fpt> u2_sinM = u2->CalcHilbertProcess(time);

    std::vector<fpt> u3_cos = u3->CalcProcess(time);

    for(int i=0;i<time.size();i++){
      re_sq = std::pow(u1_cos[i]+u2_cos[i]+p*std::cos(2*M_PI*fp*time[i]+theta_p),2);
      im_sq = std::pow(u1_sin[i]+u2_sinM[i]+p*std::sin(2*M_PI*fp*time[i]+theta_p),2);

      lgn=std::exp(u3_cos[i]*sig3+m3);

      res.push_back((std::sqrt(re_sq+im_sq))*lgn);
    }

    return res;
  }

  //fpt RiceFactor();

};


#endif
