#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "../../Models/SuzukiChannelII.hpp"

using namespace std;

int main(){

  float ti=0, tf=3, dt=0.0005;
  const short N1 = 25, N3 = 15; //number of multipath
  float sig=0.7697, fmax=91, ko=0.4045, sig3=0.0062,kc=1.735,m3=-0.3861;
  float fcut = std::sqrt(std::numbers::ln2_v<float>)*fmax;
  float p=1.567, theta_p=127*M_PI/180, theta_o=45*M_PI/41;
  float tc=1/fmax;
  

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");

  SuzukiChannelII<N1,N3,float> nt(sig,fmax,ko,theta_o,sig3,fcut,kc,m3,p,0,theta_p);
  vector<float> time((tf-ti)/dt);

  for(float i=0;i<(tf-ti)/dt;i++) time[i] = ti + i*dt;

  vector<float> ntvalue=nt.CalcChannel(time);

  OutFile.open("SuzukiChannelII.csv", ios::out | ios::trunc);

  OutFile<<"time"<<';';
  OutFile<<"nt"<<"\n";

  for(float i=0;i<(tf-ti)/dt;i++){
    OutLine.push_back(to_string(time[i]));
    OutLine.push_back(to_string(ntvalue[i]));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();

  cout<<"coherence time: "<<tc<<"\n";

  return 0;
}
