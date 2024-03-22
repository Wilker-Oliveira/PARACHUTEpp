#include <iostream>
#include "../Models/MEDSModel.hpp"

using namespace std;

int main(){

  const short N = 10000; //number of multipath
  float fc=2e9, mtv=30/3.6;
  float fmax=(mtv*fc)/3e8;
  float tc=1/fmax;
  float ti=0, tf=30*tc, dt=0.0005;

  vector<float> time((tf-ti)/dt);
  float meanPower=0;
  float meanValue=0;

  float test=0;

  MEDSModel<N,float> u1(0.7071,fmax);//0.7071 is the standard deviation of the process(sig)
  for(float i=0;i<(tf-ti)/dt;i++) time[i] = ti + i*dt;
  vector<float> u1value = u1.CalcProcess(time);

  meanPower = u1.CalcMeanPower(); 
  meanValue = u1.CalcMeanValue();


  for(float i=0;i<(tf-ti)/dt;i++) test+=u1value[i];
  test=test/((tf-ti)/dt);

  cout<<"The gaussian process u1 has the mean value of "<<meanValue<<".\n";
  cout<<"The gaussian process u1 has the empirical mean value of "<<test<<".\n";

  test=0;
  for(float i=0;i<(tf-ti)/dt;i++) test+=u1value[i]*u1value[i];
  test=test/((tf-ti)/dt);

  cout<<"The gaussian process u1 has the mean power of "<<meanPower<<".\n";
  cout<<"The gaussian process u1 has the empirical mean power of "<<test<<".\n";

  return 0;
}
