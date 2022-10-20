/** Implementation of the Process Elementary Properties class */

#ifndef PROCESSELEMETARYPROPERTIES_HPP
#define PROCESSELEMETARYPROPERTIES_HPP

#include <array>
#include <cmath>

template<short N, typename fpt>

/**
 * ProcessElementaryProperties class
 * Implements the elementary properties of the Sum of Sinusoids model, described
 * in Mobile Radio Channels (2nd edition) by Matthias Patzold, Chapter 4, section 4.2.
 */


class ProcessElementaryProperties{
    public:
    fpt CalcMeanValue();
    fpt CalcMeanPower(std::array<fpt,N> &pathGains);
    fpt CalcACF(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies, fpt tau);
    fpt CalcPSD(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies, fpt var_f);
    fpt CalcDopplerSpread(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies);
    fpt CalcPeriodicity(std::array<fpt,N> &dopplerFrequencies);//TBD
    fpt CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/);
  };

/**
 * Function for the mean value of the process
 * type: fpt
 * No parameters
 * @return 0, mean of the gaussian random process
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanValue(){return 0;}


/**
 * Function for the mean power of the process
 * type: fpt
 * @param pathGains a vector containing the path gains of the channel
 * @return MeanPower
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanPower(std::array<fpt,N> &pathGains){
    fpt MeanPower=0;
    for(int i=0; i<N; i++){
      MeanPower+=pow(pathGains[i],2)/2;
    }
    return MeanPower;
}

// finish now

/**
 * Function for the autocorrelation of the process
 * type: fpt
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @param tau, the parameter of the autocorrelation
 * @return the autocorrelation value for a determinate tau
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcACF(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies, fpt tau){
  fpt r_mu = 0;
  for(int i=0; i<N; i++){
    r_mu += pow(pathGains[i],2)*cos(2*M_PI*dopplerFrequencies[i]*tau)/2;
  }
  return r_mu;
}

/**
 * Function for the power spectral density of the process
 * type: fpt
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @param var_f, the parameter of the PSD (frequency domain)
 * @return the PSD value for a determinate var_f
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPSD(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies, fpt var_f){

  fpt S_f = 0;
  for(int i=0; i<N; i++){
    if(var_f == dopplerFrequencies[i] || var_f == -dopplerFrequencies[i]){
      S_f += pow(pathGains[i],2)/4;
    }
  }
  return S_f;
}

/**
 * Function for the Doppler spread of the process
 * type: fpt
 * @param pathGains a vector containing the path gains of the channel
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @return Doppler spread of the process
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcDopplerSpread(std::array<fpt,N> &pathGains, std::array<fpt,N> &dopplerFrequencies){

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
 * Function for the periodicity of the process
 * type: fpt
 * @param dopplerFrequencies a vector containing the Doppler frequencies of the channel
 * @return if exists, returns the periodicity. If not, returns 0.
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPeriodicity(std::array<fpt,N> &dopplerFrequencies){
}// finish later




/**
 * Function for the Probability Density of the process
 * type: vector
 * @param pathGains a vector containing the path gains of the channel
 * @return 
 */
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcPDF(fpt pathGain, fpt x /*calc the probability of u_i,n being x*/){

  if(fabs(x) < pathGain){
    return 1/(M_PI*fabs(pathGain)*pow(1-pow(x/pathGain,2),1/2) );
  }
  else return 0;
}

#endif
