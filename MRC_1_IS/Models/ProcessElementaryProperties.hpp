/** @brief Implementation of the Process Elementary Properties class. */

#ifndef PROCESSELEMETARYPROPERTIES_HPP
#define PROCESSELEMETARYPROPERTIES_HPP

#include <array>
#include <numbers>
#include <cmath>

template<short N, typename fpt>

/**
 * @brief The ProcessElementaryProperties class
 * implements the elementary properties of the Sum of Sinusoids model, described
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 4, section 4.2.
 */


class ProcessElementaryProperties{
    public:
    fpt CalcMeanValue() const;
    fpt CalcMeanPower(const std::array<fpt,N> &pathGains) const;
    fpt CalcACF(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies, fpt tau) const;
    fpt CalcPSD(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies, fpt var_f, fpt eps) const;
    fpt CalcDopplerSpread(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies) const;
    fpt CalcPeriodicity(const std::array<fpt,N> &dopplerFrequencies) const;//TBD
    fpt CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/) const;
  };

/**
 * @brief Function for the mean value of the process.
 * @brief Type: fpt.
 * @brief There are no parameters.
 * @return the value 0, theoretical mean of the gaussian random process.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanValue() const {return 0;}


/**
 * @brief Function for the mean power of the process.
 * @brief Type: fpt.
 * @param pathGains a vector containing the path gains of the channel
 * @return the theoretical mean power of the process.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanPower(const std::array<fpt,N> &pathGains) const {
    fpt MeanPower=0;
    for(int i=0; i<N; i++){
      MeanPower+=pow(pathGains[i],2)/2;
    }
    return MeanPower;
}

/**
 * @brief Function for the autocorrelation of the process.
 * @brief Type: fpt.
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @param tau the parameter of the autocorrelation
 * @return the theoretical autocorrelation value for a determinated tau.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcACF(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies, fpt tau) const {
  fpt r_mu = 0;
  for(int i=0; i<N; i++){
    r_mu += pow(pathGains[i],2)*cos(2*M_PI*dopplerFrequencies[i]*tau)/2;
  }
  return r_mu;
}

/**
 * @brief Function for the power spectral density of the process.
 * @brief Type: fpt.
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @param var_f the parameter of the PSD (frequency domain)
 * @return the theoretical PSD value for a determinate var_f.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPSD(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies, fpt var_f, fpt eps = 0.0001) const {

  fpt S_f = 0;
  for(int i=0; i<N; i++){
    if(std::fabs(var_f - dopplerFrequencies[i]) <= eps || std::fabs(dopplerFrequencies[i] + var_f) <= eps){
      S_f += pow(pathGains[i],2)/4;
    }
  }
  return S_f;
}

/**
 * @brief Function for the Doppler spread of the process.
 * @brief Type: fpt.
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @return Doppler spread of the process.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcDopplerSpread(const std::array<fpt,N> &pathGains, const std::array<fpt,N> &dopplerFrequencies) const {

  fpt dopplerSpread;
  fpt beta_i = 0, sum = 0;
  fpt sigma = sqrt(CalcMeanPower(pathGains));

  for(int i = 0; i<N; i++){
    sum += pow(pathGains[i]*dopplerFrequencies[i],2);
  }
  beta_i = 2*pow(M_PI,2)*sum;
  dopplerSpread = sqrt(beta_i)/(2*M_PI*sigma);

  return dopplerSpread;
}

/**
 * @brief Function for the periodicity of the process.
 * @brief Type: fpt.
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @return if exists, returns the periodicity. If not, returns 0.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPeriodicity(const std::array<fpt,N> &dopplerFrequencies) const {
}// finish later




/**
 * @brief Function for the Probability Density of the process.
 * @brief Type: dpt.
 * @param pathGain a fpt containing the path gains of the channel
 * @param x targeted value of the process
 * @return the PDF of the process.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/) const {

  if(fabs(x) < pathGain){
    return 1/(M_PI*fabs(pathGain)*pow(1-pow(x/pathGain,2),1/2) );
  }
  else return 0;
}

#endif
