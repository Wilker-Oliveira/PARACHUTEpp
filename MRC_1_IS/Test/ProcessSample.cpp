#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>

//Model used in this code
#include "../Models/MEAModel.hpp"

using namespace std;

int main(){

  const short N = 20; //number of multipath
  float fc=2e9, mtv=10/3.6;
  float fmax=(mtv*fc)/3e8; //for jakes PSD 
  float tc=1/fmax;
  float ti=0, tf=10*tc, dt=0.0005;
  
  float fcut = std::sqrt(std::numbers::ln2_v<float>)*fmax; //for gaussian PSD
  float kc = 2.0/std::sqrt(2.0/std::numbers::ln2_v<float>); //for gaussiain PSD

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");
  vector<float> time((tf-ti)/dt);

  MEAModel<N,float> u1(0.7071,fmax);
  MEAModel<N+1,float>  u2(0.7071,fmax);


  for(float i=0;i<(tf-ti)/dt;i++) time[i] = ti + i*dt;

  vector<float> u1value = u1.CalcProcess(time);
  vector<float> u2value = u2.CalcProcess(time);

  OutFile.open("RayleighChannelMC.csv", ios::out | ios::trunc);

  OutFile<<"time"<<';';
  OutFile<<"u1"<<';';
  OutFile<<"u2"<<"\n";


  for(float i=0;i<(tf-ti)/dt;i++){
    OutLine.push_back(to_string(time[i]));
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
