#include "SpectralProperties.hpp"

//the Oth-order Bessel function of the first kind.

//Jakes Equivalent part 

double SpectralProperties::GaussianPSD(long long freq, double std_dev, double fc){
  double var=0, AmpConst=0, exp_part=0;
  
  var=pow(std_dev,2);
  AmpConst=sqrt(log(2)/M_PI);
  exp_part=exp(-log(2)*pow(freq/fc,2));

  return (var/fc)*AmpConst*exp_part;
}
