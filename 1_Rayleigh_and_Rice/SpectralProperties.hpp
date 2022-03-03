//implementation of the Spectral properties of Rice Process's class

#ifndef SPECTRALPROPERTIES_HPP
#define SPECTRALPROPERTIES_HPP

class SpectralProperties{

public:

  //also inplement: the Oth-order Bessel function of the first kind.

  static double jakesPSD(long long int freq, double std_dev, double fmax /*Maximum Doppler frequency*/);
  static double jakesACF(short int tau, double std_dev, double fmax /*Maximum Doppler frequency*/);
  static double jakesBeta(double std_dev, double fmax /*Maximum Doppler frequency*/);

  static double GaussianPSD(long long int freq, double std_dev, double fc /*3-dB-cut-off frequency*/);
  static double GaussianACF(short int tau, double std_dev, double fc /*3-dB-cut-off frequency*/);
  static double GaussianBeta(double std_dev, double fc /*3-dB-cut-off frequency*/);
};


#endif
