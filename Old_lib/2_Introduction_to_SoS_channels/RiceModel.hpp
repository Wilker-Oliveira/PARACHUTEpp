//Implementation of the Rice model class

#ifndef RICEMODEL_HPP
#define RICEMODEL_HPP

#include <complex>

class RiceModel{

private:

  short MultPath;
  double *PathGains;
  double *DopplerFrequencies;
  double *Phases;

public:
  
  RiceModel();//set the attributes to nullptr
  ~RiceModel();//deletes the attributes

  virtual void ModelEstimate(short L/*Number of multipaths*/, float K /*Rice factor in dB*/, double fmax) = 0;

  std::complex<double> * Channel(float *t);

};

#endif
