#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <complex>
#include "../Models/MSEModel.hpp"

using namespace std;

vector<complex<float>> autoCorr(vector<complex<float>> u, float lag){
  complex<float> aux (0.0, 0.0);
  vector<complex<float>> ac;
  //vector<complex<float>> h;
  //copy(u.begin(), u.end(), back_inserter(h));
  //reverse(h.begin(), h.end());

  for(float i = 0; i < lag; i++){
    for(float j = 0; j < lag-i; j++){
      aux += u[j+i]*conj(u[j]);
    }
    ac.push_back(aux);
    aux = 0;
  }
  return ac;
}

vector<float> envACF(float sig, vector<float> acf){
  vector<float> en_acf;
  float aux;
  
  for(float i = 0; i< acf.size(); i++){
    //aux = pow(sig, 2)*(M_PI_2)*(1 + (std::pow<float>(acf[i],2)/(4*pow(sig,4))));
    aux = 2*acf[i];
    en_acf.push_back(aux);
  }
  return en_acf;
}

int main(){

  float ti=0, tf=2, dt=0.0005;
  const short N = 20, nMean=1000; //number of multipath and mean realization
  float fc=2e9, mtv=30/3.6;
  float fmax=(mtv*fc)/3e8;
  float tc=1/fmax;

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");

  MSEModel<N,float> u1(0.7071,fmax);
  MSEModel<N+1,float> u2(0.7071,fmax);
  vector<float> time((10*tc)/(dt*2));
  vector<float> T_ACFPvec((10*tc)/(dt*2));

  //Calculation of theoretical Auto correlation function
  //and generating the time vactor
  for(float i=0;i<(10*tc)/(dt*2);i++){
    time[i] = ti + i*dt;
    T_ACFPvec[i]=u2.CalcACF(time[i]);
  }

  vector<float> T_ACFvec = envACF(0.7071, T_ACFPvec);
  vector<float> u1value = u1.CalcProcess(time);
  vector<float> u2value = u2.CalcProcess(time);
  vector<complex<float>> S_ACFvec,aux;
  vector<complex<float>> uvalue;
  complex<float> ch;

  //Openning the .csv file 
  OutFile.open("RayleighChannelMSEM.csv", ios::out | ios::trunc);

  OutFile<<"time"<<';';
  OutFile<<"Theoretical_ACF"<<';';
  OutFile<<"Empiric_ACF"<<"\n";


  //Calculate the empiric Autocorrelation function for the firsh time
  for(int i = 0; i < u1value.size(); i++){
    ch = complex<float>(u1value[i], u2value[i]);
    uvalue.push_back(ch);    
  }

  S_ACFvec = autoCorr(uvalue,(10*tc)/(dt*2));
  uvalue.clear();

  //Taking the mean of the Auto correlation function

  for(int i =0;i<nMean;i++){
    u1.genPhases();
    u2.genPhases();
    u1value = u1.CalcProcess(time);
    u2value = u2.CalcProcess(time);
  
    for(int i = 0; i < u1value.size(); i++){
      ch = complex<float>(u1value[i], u2value[i]);

      uvalue.push_back(ch);
    
    }
    aux = autoCorr(uvalue,(10*tc)/(dt*2));
    uvalue.clear();
    for(float j=0;j<(10*tc)/(dt*2);j++){
      S_ACFvec[j]= (S_ACFvec[j]+aux[j]);
    }
    aux.clear();
  }

  for(float i=0;i<(10*tc)/(dt*2);i++){
    OutLine.push_back(to_string(time[i]));
    OutLine.push_back(to_string(T_ACFvec[i]));
    OutLine.push_back(to_string(real(S_ACFvec[i])/((nMean+1)*S_ACFvec.size())));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();

  cout<<"coherence time: "<<tc<<"\n";
  return 0;
}
