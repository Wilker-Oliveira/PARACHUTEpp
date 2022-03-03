//Implementation of the Rice Process class

#ifndef RICEPROCESS_HPP
#define RICEPROCESS_HPP

#include <complex>
#include "SpectralProperties.hpp"

class RiceProcess{

private:

  SpectralProperties SPtools;

public: 
  RiceProcess();
  ~RiceProcess();

  virtual std::complex<double> * RiceChannel(float *t, float K /*Rice factor in dB*/) = 0;
 
  //change char to an enum, but modify the enum to use only 1 byte instead of 4 bytes.
  double PSD(char,long long int,double,double);
  double ACF(char,short int,double,double);

  double RicePDF(double x/*the x parameter of pdf(x)*/, double std_dev, double P /*The amplitude of line-of-sight componnent*/);
  //also change this char for enum
  double LevelCrossingRate(char, double r/*Level r*/, double P/*The amplitude of line-of-sight componnent*/);
  //only for low level of r. (r<<1)
  double AverageFadingDuration(char, double r);

};

#endif
