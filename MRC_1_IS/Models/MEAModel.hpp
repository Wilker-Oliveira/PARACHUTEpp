#include "SoSModel.hpp"

template<short N, typename fpt>

class MEAModel: public SoSModel<N,fpt>{

protected: 

  //template function for gaussian PSD equation
  fpt GaussianPSDe(fpt n, fpt dfreq, float &fc){
    return ( (n/N)-std::erf(dfreq*(std::sqrt<fpt>( std::numbers::ln2_v<fpt> )/fc) ) );
  }

public:

  typedef fpt (*RTFF)(fpt,fpt,float&);

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
  fpt bissecMethod(fpt a, fpt b, fpt tol/*tolerance*/, RTFF f, short lim=100);

  //DefineModel for jakes PSD
  void DefineModel(float sig /**< std_dev in lin. */, float fmax){}

  //DefineModel for Gaussian PSD
  void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc){}


};

template<short N, typename fpt>
fpt MEAModel<N, fpt>::bissecMethod(fpt a, fpt b, fpt tol, RTFF f, short lim){};
