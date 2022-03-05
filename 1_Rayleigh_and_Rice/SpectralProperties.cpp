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

double SpectralProperties::GaussianACF(short tau, double std_dev, double fc){
  double var=0, ArgConst=0;

  var=pow(std_dev, 2);
  ArgConst=(M_PI*fc*tau)/sqrt(log(2));

  return var*exp(-pow(ArgConst, 2));
}

double SpectralProperties::GaussianBeta(double std_dev, double fc){

  return 2*pow(M_PI*fc*std_dev, 2)/log(2);
}
