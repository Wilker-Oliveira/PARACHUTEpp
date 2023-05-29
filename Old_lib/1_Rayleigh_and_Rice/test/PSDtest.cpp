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
  SpectralProperties PSD;
  
  OutFile.open("PSD_test.csv",ios::out | ios::trunc);

  OutFile<<"frequency"<<';';
  OutFile<<"JakesPSD"<<';';
  OutFile<<"GaussianPSD"<<"\n";

  //Considering standard deviation as 1.
  for(int freq=-200;freq<=200;freq++){
    OutLine.push_back(to_string(freq));
    OutLine.push_back(to_string(PSD.jakesPSD(freq, 1, fmax)));
    OutLine.push_back(to_string(PSD.GaussianPSD(freq, 1, fc)));
    copy(OutLine.begin(),OutLine.end()-1,out_it);
    OutFile<<OutLine.back()<<'\n';
    OutLine.clear();
  }
  OutFile.close();


}
