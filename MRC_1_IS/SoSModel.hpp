//Implementation of the SoS model class

#ifndef SOSMODEL_HPP
#define SOSMODEL_HPP

#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <cmath>

template<short N, typename fpt>

/**
 * SoSModel class
 * Implements the Sum of Sisoids model using the methods
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 5.
 */

/** Add a restriction for fpt, to only float and double types. */

class SoSModel{

private:

  /** short MultPath.
   * Introducing the necessary variables from equation (4.4) in Chapter 4.
   */
  std::array<fpt,N> PathGains;
  std::array<fpt,N> DopplerFrequencies;
  std::array<fpt,N> Phases;

public:
  
  SoSModel();   /**< Set the attributes to nullptr. */
  ~SoSModel();  /**< Deletes the attributes. */

  std::vector<fpt> * CalcProcess(float *t);

  /** Jakes form */
  virtual void EstimateModel(float sig /**< std_dev in lin. */, float fmax) = 0;

  /** Gaussian form */
  virtual void EstimateModel(float sig /**< std_dev in lin. */, float fc, float kc) = 0;

  void genPhases();   /**< Considering only random generation until this moment. */
  

};

template<short N, typename fpt>

/**
 * Defining the function genPhases()
 * type: void
 * No parameters
 * @return A uniformely distributed random value in the interval [0,2*pi]
 */

void SoSModel<N,fpt>::genPhases() {

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<fpt> distribution(0.0,2*M_PI);

  for(short i=0;i<N;i++) Phases[i] = distribution(generator);
}

#endif
