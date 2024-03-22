#ifndef LIBDISCRETEMODELS_H
#define LIBDISCRETEMODELS_H

//basic libs
#include <cmath>
#include <cstddef>
#include <vector>
#include <algorithm>
//random utilities
#include <random>
#include <chrono>
//Error handling
#include <stdexcept>


template <typename fpt>

/**
 * @brief The Discreteffchannels class implements the probability density functions for
 * the distributions described in the chapters 2.1 and 6.1 of  "Mobile Radio Channels" by Matthias Patzold.
 */
class DiscreteFFChannels{

private:
  std::mt19937 rdevgen; /**< @brief Define the pseudorandom numbers generator.*/
  
public:


  DiscreteFFChannels(); /**< @brief The constructor of the class Discreteffchannels. Set the attributes to nullptr.*/

/**
 * @brief Funtion that generates samples of the Gaussian distribution.
 * @param NofSamples a integer, the number of samples to be generated.
 * @param mu a fpt, the mean of the distribution.
 * @param sig a fpt, the standard deviation of the distribution.
 * @return a vector of length Nofsamples with the samples of the Gaussian distribution.  
 */
  std::vector<fpt> gaussDistribution(int NofSamples, fpt mu, fpt sig);

/**
 * @brief Funtion that generates samples of the Lognormal distribution.
 * @param NofSamples a integer, the number of samples to be generated.
 * @param mu a fpt, the mean of the distribution.
 * @param sig a fpt, the standard deviation of the distribution.
 * @return a vector of length Nofsamples with the samples of the Lognormal distribution.
 */
  std::vector<fpt> lognormDistribution(int NofSamples, fpt mu, fpt sig);

/**
 * @brief Funtion that generates samples of the Rayleigh distribution.
 * @param NofSamples a integer, the number of samples to be generated.
 * @param sig a fpt, the standard deviation of the distribution.
 * @return a vector of length Nofsamples with the samples of the Rayleigh distribution.
 */
  std::vector<fpt> rayleighDistribution(int NofSamples, fpt sig);

/**
 * @brief Funtion that generates samples of the Rice distribution.
 * @param NofSamples a integer, the number of samples to be generated.
 * @param sig a fpt, the standard deviation of the distribution.
 * @param rho a fpt, the amplitude of the line-of-sight channel component.
 * @return a vector of length Nofsamples with the samples of the Rice distribution.
 */
  std::vector<fpt> riceDistribution(int NofSamples, fpt sig, fpt rho);

/**
 * @brief Funtion that generates samples of the Suzuki distribution.
 * @param NofSamples a integer, the number of samples to be generated.
 * @param sig0 a fpt, the standard deviation of the Rice distribution.
 * @param rho a fpt, the amplitude of the line-of-sight channel component.
 * @param mu a fpt, the mean of the Lognormal distribution.
 * @param sigu a fpt, the standard deviation of the Lognormal distribution.
 * @return a vector of length Nofsamples with the samples of the Suzuki distribution.
 */  
  std::vector<fpt> suzukiDistribution(int NofSamples, fpt sig0, fpt rho, fpt mu, fpt sigu);
  
};

template <typename fpt>
/**
 * @details
 * Performs the restriction for the fpt template parameter, must be only
a Floating Point type.
 * @brief Sets the seed for the pseudorandom numbers. 
 */
DiscreteFFChannels<fpt>::DiscreteFFChannels(){
  if(!std::is_floating_point<fpt>()) throw std::runtime_error("fpt should be a Floating Point Type");
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  this->rdevgen.seed(seed);
}

template <typename fpt>
/**
 * @details
 * The gaussdistribution() function generates random samples that follow
a Gaussian distribution, described by the mean and the standard deviation.
 * @details
 * The PDF of the Gaussian distribution is described by equation (2.30) in "Mobile
Radio Channels" by Matthias Patzold.
 */
std::vector<fpt> DiscreteFFChannels<fpt>::gaussDistribution(int NofSamples, fpt mu, fpt sig){
  std::vector<fpt> ret(NofSamples,0);
  std::normal_distribution<fpt> distribution(mu,sig);
  
  std::generate(ret.begin(), ret.end(), [&]() mutable {return distribution(rdevgen);});

  return ret;
}

template <typename fpt>
/**
 * @details
 * The lognormdistribution() function generates random samples that follow
a Lognormal distribution, described by the mean and the standard deviation in a
logarithmic base.
 * @details
 * The PDF of the Lognormal distribution is described by equation (2.51) in "Mobile
Radio Channels" by Matthias Patzold.
 */
std::vector<fpt> DiscreteFFChannels<fpt>::lognormDistribution(int NofSamples, fpt mu, fpt sig) {
  std::vector<fpt> ret(NofSamples,0);
  std::lognormal_distribution<double> distribution(mu, sig);

  std::generate(ret.begin(), ret.end(), [&]() mutable {return distribution(rdevgen);});

  return ret;
}

template <typename fpt>
/**
 * @details
 * The rayleighdistribution() function generates random samples that follow
a Rayleigh distribution, described by the standard deviation.
 * @details
 * The PDF of the Rayleigh distribution is described by equation (2.40) in "Mobile
Radio Channels" by Matthias Patzold. This implementation uses the Rayleigh distributed random number sequence function in example 7-20 of "Probability, Random Variables and Stochastic Processes" by Athanasios Papoulis, in which the Rayleigh distributed sequence is generated using the Percentile Transformation Method, with the natural logarithm of a uniform random variable.
 */
std::vector<fpt> DiscreteFFChannels<fpt>::rayleighDistribution(int NofSamples, fpt sig){
  std::vector<fpt> ret(NofSamples,0);
  std::uniform_real_distribution<fpt> distribution(0.0, 1.0);

  std::generate(ret.begin(), ret.end(), [&]() mutable { return sig*std::pow(-2*std::log(distribution(rdevgen)),0.5);});

  return ret;
}

template <typename fpt>
/**
 * @details
 * The ricedistribution() function generates random samples that follow
a Rice distribution, described by the standard deviation and the amplitude of
the line-of-sight channel component.
 * @details
 * The PDF of the Rice distribution is described by equation (2.44)
in "Mobile Radio Channels" by Matthias Patzold.
 */
std::vector<fpt> DiscreteFFChannels<fpt>::riceDistribution(int NofSamples, fpt sig, fpt rho){
  std::vector<fpt> ret(NofSamples,0);  
  std::normal_distribution<fpt> distribution1(0.0,sig);

  std::generate(ret.begin(), ret.end(), [&]() mutable { return std::sqrt(std::pow((distribution1(rdevgen)+rho),2)+std::pow(distribution1(rdevgen),2));} );
  
  return ret;
}

template <typename fpt>
/**
 * @details
 * The suzukidistribution() function generates random samples that follow
a Suzuki distribution, described by the product of a Rice distributed variable
and a Lognormal distributed variable, as shown in equation (6.54) of "Mobile
Radio Channels" by Matthias Patzold. The Suzuki distribution is described by the
standard deviation of the Rice distribution, the amplitude of the line-of-sight
channel component, and the mean and standard deviation of the lognormal
distribution.
 * @details
 * The PDF of the Suzuki distribution is described by equation (6.56) in "Mobile Radio Channels" by Matthias Patzold.
 */
std::vector<fpt> DiscreteFFChannels<fpt>::suzukiDistribution(int NofSamples, fpt sig0, fpt rho, fpt mu, fpt sigu){
  std::vector<fpt> ret;
  ret.reserve(NofSamples);

  std::vector<fpt> lamb = lognormDistribution(NofSamples,mu,sigu);
  std::vector<fpt> zeta = riceDistribution(NofSamples,sig0,rho);

  for (std::size_t i = 0; i < NofSamples; ++i)
    ret[i] = lamb[i]*zeta[i];
  
  return ret;
  }



#endif
