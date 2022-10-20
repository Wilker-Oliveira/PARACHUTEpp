#include <iostream>
#include "../Models/MEDModel.hpp"

using namespace std;

int main(){

  float ti=0, tf=2, dt=0.0005;
  const short N = 5; //number of multipath
  float fc=2e9, mtv=30/3.6;
  float fmax=(mtv*fc)/3e8;
  float tc=1/fmax;

  float meanPower=0;
  float meanValue=0;

  MEDModel<N,float> u1;
  u1.genPhases();
  u1.DefineModel(0.7071,fmax);//0.7071 is the standard deviation of the process(sig)

  meanPower = u1.CalcMeanPower(); 
  meanValue = u1.CalcMeanValue();

  cout<<"The process gaussian process u1 has the mean power of "<<meanPower<<" and the mean value of "<<meanValue<<".\n";
  
  //mean power is ideally equal to sig^2

  return 0;
}
