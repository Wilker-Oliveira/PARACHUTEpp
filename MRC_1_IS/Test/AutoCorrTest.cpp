#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <complex>

/** Model used in this code.*/
#include "../Models/JakesModel.hpp"

using namespace std;

/**
 * @brief This function computes the empirical autocorrelation of the process \f$\mu\f$.
 * @param u a vector of float complex numbers, the simulated process
 * @param lag a float, the maximum \f$\tau\f$ parameter of the autocorrelation
 * @return a vector containing the autocorrelation of \f$\mu\f$ for each time difference in \f$[0,\tau]\f$
 */
vector<complex<float>> autoCorr(vector<complex<float>> u, float lag){
  complex<float> aux (0.0, 0.0);
  vector<complex<float>> ac;

  for(float i = 0; i < lag; i++){
    for(float j = 0; j < lag-i; j++){
      aux += u[j+i]*conj(u[j]);
    }
    ac.push_back(aux);
    aux = 0;
  }
  return ac;
}

/**
 * @brief This function computes the theoretical autocorrelation of the complex envelope \f$\mu\f$ accordingly to the equation (3.24) in chapter 3 of Mobile Radio Channels by Matthias Patzold.
 * @param sig a float, the standard deviation of the process.
 * @param acf a vector of floats, the theoretical autocorrelation of the process \f$\mu_i\f$
 * @return a vector containing the autocorrelation of \f$\mu\f$ for each time difference in \f$[0,\tau]\f$.
 */
vector<float> envACF(float sig, vector<float> acf1, vector<float> acf2){
  vector<float> en_acf;
  float aux;
  
  for(float i = 0; i< acf1.size(); i++){
    aux = acf1[i] + acf2[i];
    en_acf.push_back(aux);
  }
  return en_acf;
}

int main(){

  /** Declaring the variables */
  const short N = 20, nMean=1000; //number of multipath and mean realization
  float fc=2e9, mtv=30/3.6;
  float fmax=(mtv*fc)/3e8; // for Jakes PSD
  float tc=1/fmax;
  float ti=0, tf=2, dt=0.0005;

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");

  /** Time vectors */
  vector<float> time((10*tc)/(dt*2));       // for the process simulation
  vector<float> T_ACFPvec1((10*tc)/(dt*2));  // for the autocorrelation
  vector<float> T_ACFPvec2((10*tc)/(dt*2));  // for the autocorrelation

  /** Processes */
  JakesModel<N,float> u1(0.7071,fmax,false);
  JakesModel<N+1,float>  u2(0.7071,fmax,true);

  /** Calculation of theoretical autocorrelation function of the process \f$\mu_i\f$ and generating the time vector. */
  for(float i=0;i<(10*tc)/(dt*2);i++){
    time[i] = ti + i*dt;
    T_ACFPvec1[i]=u1.CalcACF(time[i]);
    T_ACFPvec2[i]=u2.CalcACF(time[i]);
  }

  /** Empirical autocorrelation of the complex envelope \f$\mu\f$*/
  vector<float> T_ACFvec = envACF(0.7071, T_ACFPvec1, T_ACFPvec2);

  vector<float> u1value = u1.CalcProcess(time);
  vector<float> u2value = u2.CalcProcess(time);

  vector<complex<float>> S_ACFvec;    // vector to store the empirical autocorrelation values
  vector<complex<float>> uvalue,aux; 
  complex<float> ch;

  //Calculate the empiric Autocorrelation function for the firsh time
  for(int i = 0; i < u1value.size(); i++){
    ch = complex<float>(u1value[i], u2value[i]);
    uvalue.push_back(ch);    
  }

  S_ACFvec = autoCorr(uvalue,(10*tc)/(dt*2));
  uvalue.clear();

  /**Taking the mean of the empirical autocorrelation function.*/

  for(int i =0;i<nMean;i++){
    // Generate the random phases
    u1.genPhases();
    u2.genPhases();

    // Compute the processes
    u1value = u1.CalcProcess(time);
    u2value = u2.CalcProcess(time);

    // Generate the complex envelope
    for(int i = 0; i < u1value.size(); i++){
      ch = complex<float>(u1value[i], u2value[i]);
      uvalue.push_back(ch);
    }
    
    // Compute the empirical autocorrelation in auxiliar vector
    aux = autoCorr(uvalue,(10*tc)/(dt*2));

    // Clear complex envelope vector
    uvalue.clear();

    // Store the autocorrelation function result
    for(float j=0;j<(10*tc)/(dt*2);j++){
      S_ACFvec[j]= (S_ACFvec[j]+aux[j]);
    }

    // Clear auxiliar vector
    aux.clear();
  }

  //Openning the .csv file 
  OutFile.open("RayleighChannelJakes.csv", ios::out | ios::trunc);

  OutFile<<"time"<<';';
  OutFile<<"Theoretical_ACF"<<';';
  OutFile<<"Empirical_ACF"<<"\n";

  for(float i=0;i<(10*tc)/(dt*2);i++){
    OutLine.push_back(to_string(time[i]));
    OutLine.push_back(to_string(T_ACFvec[i]));
    OutLine.push_back(to_string(real(S_ACFvec[i])/((nMean+1)*S_ACFvec.size())));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();

  
  return 0;
}
