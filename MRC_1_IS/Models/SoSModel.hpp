/**
 * @brief Implementation of the SoS model class.
*/

#ifndef SOSMODEL_HPP
#define SOSMODEL_HPP

#include <vector>
#include <random>
#include <chrono>
#include <stdexcept>
#include "ProcessElementaryProperties.hpp"

template<short N, typename fpt>

/**
 * @brief The SoSModel class
 * implements the Sum of Sinusoids model using the methods
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 5.
 */

/** @brief Restriction for fpt, to only float, double and long double types. */

class SoSModel{
  /**
 *@brief Short MultPath.
 * Introducing the necessary variables from equation (4.4) in Chapter 4.
 */

protected:

  std::array<fpt,N> pathGains; /**< N-dimensional array of type fpt containing the path gains \f$c_{i,n}\f$ of the process. */
  std::array<fpt,N> dopplerFrequencies; /**< N-dimensional array of type fpt containing the Doppler frequencies \f$f_{i,n}\f$ of the process. */
  std::array<fpt,N> phases; /**< N-dimensional array of type fpt containing the phases \f$\theta_{i,n}\f$ of the process. */

public:
  
  SoSModel();   /**< Set the attributes to nullptr. */  
  //~SoSModel();  /**< Deletes the attributes. */

  ProcessElementaryProperties<N,fpt> processProperties; /**< An object of type ProcessElementaryProperties to be used for this class and all derived classes to execute the calculation of the gaussian process referred in the class. */


  /** @brief Function implementing equation (4.4) from Mobile Radio Channels by Matthias Patzold. */
  std::vector<fpt> CalcProcess(const std::vector<float> &time) const;

  /** @brief Function that calculates the process after passing through a hilbert transform. */
  std::vector<fpt> CalcHilbertProcess(const std::vector<float> &time) const;

  /** @brief Jakes power spectral density form */
  virtual void DefineModel(float sig, float fmax) = 0;

  /** @brief Gaussian power spectral density form */
  virtual void DefineModel(float sig, float fc, float kc) = 0;

  void genPhases();   /**< @brief Considering only random generation until this moment. */
  
  /** @brief Calculating the parametric properties of the SoSModel through the processProperties class */
  fpt CalcMeanValue() const;
  fpt CalcMeanPower() const;
  fpt CalcACF(fpt tau) const;
  fpt CalcPSD(fpt f) const;
  std::vector<fpt> CalcDiscreteDopplerSpectrum() const;
  fpt CalcDopplerSpread() const;
  //add periocity after implementation
};

template<short N, typename fpt>
SoSModel<N,fpt>::SoSModel(){
  if(!std::is_floating_point<fpt>()) throw std::runtime_error("fpt should be a Floating Point Type");
  if(N<0) throw std::runtime_error("The number of multipath must be positive or zero");
}

template<short N, typename fpt>

/**
 * @brief Defining the function genPhases().
 * @brief This function generates a uniformely distributed pseudorandom value in the interval \f$[0,2*\pi]\f$.
 * @brief Type: void.
 * @brief There are no parameters.
 */

void SoSModel<N,fpt>::genPhases() {

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<fpt> distribution(0.0,2*M_PI);

  for(short i=0;i<N;i++) phases[i] = distribution(generator);
}

template<short N, typename fpt>

/**
 * @brief Defining the function CalcProcess().
 * @brief Type: vector.
 * @param time a vector representing the considered time interval.
 * @return A vector containing the channel realization for each instant of time.
 */

//return this by reference!
std::vector<fpt> SoSModel<N,fpt>::CalcProcess(const std::vector<float> &time) const {

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

template<short N, typename fpt>

/**
 * @brief Defining the function CalcHilbertProcess().
 * @brief Type: vector.
 * @param time a vector representing the considered time interval.
 * @return A vector containing the channel realization for each instant of time.
 */

//return this by reference!
std::vector<fpt> SoSModel<N,fpt>::CalcHilbertProcess(const std::vector<float> &time) const {

  fpt aux=0;
  std::vector<fpt> res;
  
  for (auto& t : time){
    for(short i=0;i<N;i++){
      aux+=pathGains[i]*sin(2*M_PI*dopplerFrequencies[i]*t + phases[i]);
    }
    res.push_back(aux);
    aux=0;
  }
  return res;
}

/**
 *@brief Implementing the function CalcMeanValue().
 *@brief Type: fpt.
 *@brief There are no parameters.
 *@return The mean value of the process.
 */
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcMeanValue() const {
  return processProperties.CalcMeanValue();
}

/**
 *@brief Implementing the function CalcMeanPower().
 *@brief Type: fpt.
 *@brief There are no parameters.
 *@return The theoretical mean power of the process.
 */
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcMeanPower() const {
  return processProperties.CalcMeanPower(this->pathGains);
}

/**
 *@brief Implementing the function CalcACF().
 *@brief Type: fpt.
 *@brief There are no parameters.
 *@return The theoretical autocorrelation function of the process.
 */
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcACF(fpt tau) const {
  return processProperties.CalcACF(this->pathGains, this->dopplerFrequencies, tau);
}

/**
 *@brief Implementing the function CalcPSD().
 *@brief Type: fpt.
 *@brief There are no parameters.
 *@return the theoretical power spectral density of the process.
 */
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcPSD(fpt f) const {
  return processProperties.CalcPSD(this->pathGains, this->dopplerFrequencies, f);
}

template<short N, typename fpt>
std::vector<fpt> SoSModel<N,fpt>::CalcDiscreteDopplerSpectrum() const {

    return processProperties.CalcDiscreteDopplerSpectrum(this->pathGains, this->dopplerFrequencies);
}


/**
 *@brief Implementing the function CalcDopplerspread().
 *@brief Type: fpt.
 *@brief There are no parameters.
 *@return The theoretical doppler spread of the process.
 */
template<short N, typename fpt>
fpt SoSModel<N,fpt>::CalcDopplerSpread() const {
  return processProperties.CalcDopplerSpread(this->pathGains, this->dopplerFrequencies);
}


#endif
