#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "../../../src/Models/SuzukiChannelI.hpp"

using namespace std;

int main(){

  float ti=0, tf=3, dt=0.0005;
  const short N1 = 20,N2 = 15, N3 = 15; //number of multipath
  float sig=0.4497, fmax=91, ko=5.9e-8, sig3=0.0101,kc=20,m3=0.0875;
  float fcut = std::sqrt(std::numbers::ln2_v<float>)*fmax;
  float p=0.9856, fp=0.7326*fmax, theta_p=0;
  float tc=1/fmax;
  

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");

  SuzukiChannelI<N1,N2,N3,float> nt(sig,fmax,ko,sig3,fcut,kc,m3,p,fp,theta_p);
  vector<float> time((tf-ti)/dt);

  for(float i=0;i<(tf-ti)/dt;i++) time[i] = ti + i*dt;

  vector<float> ntvalue=nt.CalcChannel(time);

  OutFile.open("SuzukiChannelI.csv", ios::out | ios::trunc);

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
