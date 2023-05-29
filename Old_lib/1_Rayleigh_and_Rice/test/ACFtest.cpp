#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "../SpectralProperties.hpp"

using namespace std;

int main(){

  double fmax=0,fc=0;
  fmax=91;//Hz
  fc=sqrt(log(2))*fmax;//Hz

  vector<string> OutLine;
  fstream OutFile;
  ostream_iterator<string> out_it (OutFile,";");
  SpectralProperties ACL;
  
  OutFile.open("ACF_test.csv",ios::out | ios::trunc);

  OutFile<<"tau"<<';';
  OutFile<<"JakesACF"<<';';
  OutFile<<"GaussianACF"<<"\n";

  //Considering standard deviation as 1.
  for(float tau=0;tau<=0.05;tau+=0.0005){
    OutLine.push_back(to_string(tau));
    OutLine.push_back(to_string(ACL.jakesACF(tau, 1, fmax)));
    OutLine.push_back(to_string(ACL.GaussianACF(tau, 1, fc)));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();


}
