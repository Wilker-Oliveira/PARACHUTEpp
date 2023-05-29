//implementation of the Spectral properties of Rice Process's class

#ifndef SPECTRALPROPERTIES_HPP
#define SPECTRALPROPERTIES_HPP

#include <cmath>

class SpectralProperties{

public:

  //the Oth-order Bessel function of the first kind.
  static double besselJ0(double x, float step);
  static double besselJ0(double x, int t /*Number of terms*/);

  static double jakesPSD(double freq, double std_dev, double fmax /*Maximum Doppler frequency*/);
  static double jakesACF(float tau, double std_dev, double fmax /*Maximum Doppler frequency*/);
  static double jakesBeta(double std_dev, double fmax /*Maximum Doppler frequency*/);

  static double GaussianPSD(double freq, double std_dev, double fc /*3-dB-cut-off frequency*/);
  static double GaussianACF(float tau, double std_dev, double fc /*3-dB-cut-off frequency*/);
  static double GaussianBeta(double std_dev, double fc /*3-dB-cut-off frequency*/);
};


#endif
