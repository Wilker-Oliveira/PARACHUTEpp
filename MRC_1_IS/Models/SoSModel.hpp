/** Implementation of the SoS model class */

#ifndef SOSMODEL_HPP
#define SOSMODEL_HPP

#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <stdexcept>
#include "ProcessElementaryProperties.hpp"

template<short N, typename fpt>

/**
 * SoSModel class
 * Implements the Sum of Sisoids model using the methods
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 5.
 */

/** restriction for fpt, to only float, double and long double types. */

class SoSModel{

protected:

  /** short MultPath.
   * Introducing the necessary variables from equation (4.4) in Chapter 4.
   */
  std::array<fpt,N> pathGains;
  std::array<fpt,N> dopplerFrequencies;
  std::array<fpt,N> phases;

public:
  
  SoSModel();   /**< Set the attributes to nullptr. */
  //~SoSModel();  /**< Deletes the attributes. */

  ProcessElementaryProperties<N,fpt> processProperties;


  /** Function implementing equation (4.4). */
  std::vector<fpt> CalcProcess(std::vector<float> &time);

  /** Jakes form */
  virtual void DefineModel(float sig /**< std_dev in lin. */, float fmax) = 0;

  /** Gaussian form */
  virtual void DefineModel(float sig /**< std_dev in lin. */, float fc, float kc) = 0;

  void genPhases();   /**< Considering only random generation until this moment. */
  
  /** Calculating the properties of the SoSModel through the processPRoperties class */
  fpt CalcMeanValue();
  fpt CalcMeanPower();
  fpt CalcACF(fpt tau);
  fpt CalcPSD(fpt f);
  fpt CalcDopplerSpread();
  //add periocity after implementation
};
//create another constructor to initiate all the class members at once
template<short N, typename fpt>
SoSModel<N,fpt>::SoSModel(){
  if(!std::is_floating_point<fpt>()) throw std::runtime_error("fpt should be a Floating Point Type");
  if(N<=0) throw std::runtime_error("The number of multipath must be positive");
}

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

  for(short i=0;i<N;i++) phases[i] = distribution(generator);
}

template<short N, typename fpt>

/**
 * Defining the function CalcProcess()
 * type: vector
 * @param time a vector representing the considered time interval
 * @return res a vector containing the channel realization for each instant of time
 */


//return this by reference!
std::vector<fpt> SoSModel<N,fpt>::CalcProcess(std::vector<float> &time){

  fpt aux=0;
  std::vector<fpt> res;
  
  for (auto& t : time){
    for(short i=0;i<N;i++){
      aux+=pathGains[i]*cos(2*M_PI*dopplerFrequencies[i]*t + phases[i]);
    }
    res.push_back(aux);
    aux=0;
  }
  return res;
}

//Document me
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcMeanValue(){
  return processProperties.CalcMeanValue();
}

//Document me
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcMeanPower(){
  return processProperties.CalcMeanPower(this->pathGains);
}

//Document me
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcACF(fpt tau){
  return processProperties.CalcACF(this->pathGains, this->dopplerFrequencies, tau);
}

//Document me
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcPSD(fpt f){
  return processProperties.CalcPSD(this->pathGains, this->dopplerFrequencies, f);
}

//Document me
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcDopplerSpread(){
  return processProperties.CalcDopplerSpread(this->pathGains, this->dopplerFrequencies);
}


#endif
