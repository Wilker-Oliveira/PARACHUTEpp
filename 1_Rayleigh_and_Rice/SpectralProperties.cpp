#include "SpectralProperties.hpp"

//the Oth-order Bessel function of the first kind.

//Numerical Integration form (one third of simpson)
double SpectralProperties::besselJ0(double x, float step=0.01){
  double I=0,s1=0,s2=0,s3=0;

  auto f=[](double x,double o){ return cos(x*cos(o)); };

  for(float i=step;i<(M_PI_2-step);i+=step){
    if(int(i/step)%2!=0) s2+=f(x,i);
    else s3+=f(x,i);
  }
  
  s1=f(x,0)+f(x,M_PI_2);

  I=(s1+4*s2+2*s3)*(step/3);
  
  return M_2_PI*I;
};

//Truncated infinite power-series form
double SpectralProperties::besselJ0(double x, int t=4){
  double J0=1,x2=x*x,prod=1;
  
  for(int i=1;i<=t;i++){
    prod*=pow((2*i),2);
    J0+=pow(-x2,i)/prod;
  }

  return J0;
}



//Jakes Equivalent part 
double SpectralProperties::jakesPSD(double freq, double std_dev, double fmax){
  double var=0;
  var=pow(std_dev,2);

  if(fabs(freq)<=fmax){
    if(freq==fmax)
      freq-=0.3;
    else if (-freq==fmax)
      freq+=0.3;
    return var/(M_PI*fmax*sqrt(1-pow(freq/fmax,2)));
  }
  return 0;
}

double SpectralProperties::jakesACF(float tau, double std_dev, double fmax){
 double var=0,arg=0;
 var=pow(std_dev,2);

 arg=2*M_PI*fmax*tau;

 //limtis's choose based on empirical test
 if(arg>=-3 && arg <3.1) 
   return var*besselJ0(arg,4);
 
 return var*besselJ0(arg,0.01f);

}

double SpectralProperties::jakesBeta(double std_dev, double fmax){

  return 2*pow(M_PI*fmax*std_dev, 2);
}


//Gaussian Equivalent part
double SpectralProperties::GaussianPSD(double freq, double std_dev, double fc){
  double var=0, AmpConst=0, exp_part=0;
  
  var=pow(std_dev,2);
  AmpConst=sqrt(log(2)/M_PI);
  exp_part=exp(-log(2)*pow(freq/fc,2));

  return (var/fc)*AmpConst*exp_part;
}

double SpectralProperties::GaussianACF(float tau, double std_dev, double fc){
  double var=0, ArgConst=0;

  var=pow(std_dev, 2);
  ArgConst=(M_PI*fc*tau)/sqrt(log(2));

  return var*exp(-pow(ArgConst, 2));
}

double SpectralProperties::GaussianBeta(double std_dev, double fc){

  return 2*pow(M_PI*fc*std_dev, 2)/log(2);
}
