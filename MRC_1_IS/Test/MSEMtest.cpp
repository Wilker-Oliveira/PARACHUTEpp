#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "../Models/MSEModel.hpp"

using namespace std;

//double AutoCorrelation(); 

int main(){

  float ti=0, tf=2, dt=0.0005;
  const short N = 20; //number of multipath
  float fc=2e9, mtv=30/3.6;
  float fmax=(mtv*fc)/3e8;
  float tc=1/fmax;
  

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");

  MSEModel<N,float> u1(0.7071,fmax),u2(0.7071,fmax);
  vector<float> time((tf-ti)/dt);
  vector<float> ACFvec((tf-ti)/dt);

    for(float i=0;i<(tf-ti)/dt;i++){
      time[i] = ti + i*dt;
      ACFvec[i]=u2.CalcACF(time[i]);
    }
  
  vector<float> u1value = u1.CalcProcess(time);
  vector<float> u2value = u2.CalcProcess(time);
  vector<float> aux1,aux2;

  for(int j=0; j<10000; j++){
    u1.genPhases();
    u2.genPhases();
    aux1=u1.CalcProcess(time);
    aux2=u2.CalcProcess(time);
    for(float k=0;k<(tf-ti)/dt;k++){
      u1value[k] = aux1[k] + u1value[k];
      u2value[k] = aux2[k] + u2value[k];
    }
    aux1.clear();
    aux2.clear();
  }

  for(float k=0;k<(tf-ti)/dt;k++){
    u1value[k] /= 10000;
    u2value[k] /= 10000;
  }

  
  OutFile.open("RayleighChannelMSEM.csv", ios::out | ios::trunc);

  OutFile<<"time"<<';';
  OutFile<<"ACF"<<';';
  OutFile<<"u1"<<';';
  OutFile<<"u2"<<"\n";


  for(float i=0;i<(tf-ti)/dt;i++){
    OutLine.push_back(to_string(time[i]));
    OutLine.push_back(to_string(ACFvec[i]));
    OutLine.push_back(to_string(u1value[i]));
    OutLine.push_back(to_string(u2value[i]));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();

  cout<<"coherence time: "<<tc<<"\n";
  return 0;
}
