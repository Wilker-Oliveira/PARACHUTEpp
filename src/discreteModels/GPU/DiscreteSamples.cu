#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "DiscreteModels_cuda.hpp"

int main(){

  //Number of samples
  const int N = 1000000;

  //Gaussian parameters
  float gauss_mu=0.0, gauss_sig0=std::sqrt(1.0),gauss_sig1=std::sqrt(2.0), gauss_sig2=std::sqrt(3.0);

  //Rayleigh parameters
  float ray_sig0=std::sqrt(1.0), ray_sig1=std::sqrt(2.0), ray_sig2=std::sqrt(3.0);

  //Rice parameters
  float rice_rho0 = 0.0, rice_rho1 = 1.0, rice_rho2 = 2.0, rice_sig0 = 1;

  //Lognormal parameters
  float log_sig0 = std::sqrt(1.0), log_sig1 = std::sqrt(2.0), log_sig2 = std::sqrt(3.0), log_mu = 0.0;

  //Suzuki parameters
  float suzuki_rho0 = 0.0, suzuki_rho1 = 3.0; // rho values for figures 2.4 and 6.7.
  float suzuki_sig00 = std::sqrt(1.0), suzuki_sig01 = std::sqrt(3.0); //Std values for rayleigh component of figure 2.4.
  float suzuki_sigmu00 = 0.0, suzuki_sigmu01 = std::sqrt(1.0 / 3.0), suzuki_sigmu02 = std::sqrt(1.0 / 2.0); // Std values for lognormal component of figure 2.4.
  float suzuki_mu0 = 0.0; // Mean value for lognormal component of figure 2.4.
  float suzuki_sigmu10 = 1.0 / 3.0, suzuki_sigmu11 = 2.0 / 3.0, suzuki_sigmu12 = 1.0; // Std values for lognormal component of figure 6.7.
  float suzuki_mu10 = -std::pow(suzuki_sigmu10,2)/2.0, suzuki_mu11=-std::pow(suzuki_sigmu11,2)/2.0, suzuki_mu12=-std::pow(suzuki_sigmu12,2)/2.0; // Mean values for lognormal component of figure 6.7.
  
  //auxiliary vectors to store the values
  std::vector<float> aux0, aux1, aux2, aux3, aux4, aux5;
  aux0.reserve(N);
  aux1.reserve(N);
  aux2.reserve(N);
  aux3.reserve(N);
  aux4.reserve(N);
  aux5.reserve(N);  

  //Output objects
  std::vector<std::string> OutLine;
  std::fstream OutFile;
  std::ostream_iterator<std::string> out_it(OutFile, ";");
  
  //creating the channel object
  auto channel = new gpu_model::DiscreteFFChannels<float>;

  // ================================================================ Gaussian Distribution ================================================================

  //Generating Gaussian Samples
  aux0 = channel->gaussDistribution(N,gauss_mu,gauss_sig0);
  aux1 = channel->gaussDistribution(N,gauss_mu,gauss_sig1);
  aux2 = channel->gaussDistribution(N,gauss_mu,gauss_sig2);

  //Saving Gaussian Samples to a file
  OutFile.open("GaussianSamples.csv", std::ios::out | std::ios::trunc);
  OutFile << "sig0"<< ';' ;
  OutFile << "sig1"<< ';' ;
  OutFile << "sig2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();

  // ================================================================ Rayleigh Distribution ================================================================
  
  //Generating Rayleigh Samples
  aux0 = channel->rayleighDistribution(N,ray_sig0);
  aux1 = channel->rayleighDistribution(N,ray_sig1);
  aux2 = channel->rayleighDistribution(N,ray_sig2);

    
  //Saving Rayleigh Samples to a file
  OutFile.open("RayleighSamples.csv", std::ios::out | std::ios::trunc);
  OutFile << "sig0"<< ';' ;
  OutFile << "sig1"<< ';' ;
  OutFile << "sig2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));    
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();

  // ================================================================ Rice Distribution ================================================================

  //Saving Rice Samples to a file
  aux0 = channel->riceDistribution(N,rice_sig0,rice_rho0);
  aux1 = channel->riceDistribution(N,rice_sig0,rice_rho1);
  aux2 = channel->riceDistribution(N,rice_sig0,rice_rho2);

    //Saving Rice Samples to a file
  OutFile.open("RiceSamples.csv", std::ios::out | std::ios::trunc);
  OutFile << "rho0"<< ';' ;
  OutFile << "rho1"<< ';' ;
  OutFile << "rho2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();

  // ================================================================ Lognormal Distribution ================================================================

  //Saving Lognormal Samples to a file
  aux0 = channel->lognormDistribution(N,log_mu,log_sig0);
  aux1 = channel->lognormDistribution(N,log_mu,log_sig1);
  aux2 = channel->lognormDistribution(N,log_mu,log_sig2);

  //Saving Lognormal Samples to a file
  OutFile.open("LognormalSamples.csv", std::ios::out | std::ios::trunc);
  OutFile << "sig0"<< ';' ;
  OutFile << "sig1"<< ';' ;
  OutFile << "sig2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();

  // ================================================================ Suzuki Distribution ================================================================

  // ************************************ FIGURE 2.4 ************************************
  //Saving Suzuki Samples to a file
  // sig0^2 = 1.0
  aux0 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho0,suzuki_mu0,suzuki_sigmu00);
  aux1 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho0,suzuki_mu0,suzuki_sigmu01);
  aux2 = channel->suzukiDistribution(N,suzuki_sig00, suzuki_rho0,suzuki_mu0, suzuki_sigmu02);
  // sig0^2 = 3.0
  aux3 = channel->suzukiDistribution(N,suzuki_sig01,suzuki_rho0,suzuki_mu0,suzuki_sigmu00);
  aux4 = channel->suzukiDistribution(N,suzuki_sig01,suzuki_rho0,suzuki_mu0,suzuki_sigmu01);
  aux5 = channel->suzukiDistribution(N,suzuki_sig01,suzuki_rho0,suzuki_mu0,suzuki_sigmu02);

  //Saving Suzuki Samples to a file
  OutFile.open("SuzukiSamples_figure2_4.csv", std::ios::out | std::ios::trunc);
  OutFile << "sig0_sigmu0"<< ';' ;
  OutFile << "sig0_sigmu1"<< ';' ;
  OutFile << "sig0_sigmu2"<< ';' ;
  OutFile << "sig1_sigmu0"<< ';' ;
  OutFile << "sig1_sigmu1"<< ';' ;
  OutFile << "sig1_sigmu2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));
    OutLine.push_back(std::to_string(aux3[i]));
    OutLine.push_back(std::to_string(aux4[i]));
    OutLine.push_back(std::to_string(aux5[i]));
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();
  aux3.clear();
  aux4.clear();
  aux5.clear();

  // ************************************ FIGURE 6.7 ************************************
  //Saving Suzuki Samples to a file
  // rho = 0.0
  aux0 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho0,suzuki_mu10,suzuki_sigmu10);
  aux1 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho0,suzuki_mu11,suzuki_sigmu11);
  aux2 = channel->suzukiDistribution(N, suzuki_sig00, suzuki_rho0,suzuki_mu12, suzuki_sigmu12);
  // rho = 3.0
  aux3 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho1,suzuki_mu10,suzuki_sigmu10);
  aux4 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho1,suzuki_mu11,suzuki_sigmu11);
  aux5 = channel->suzukiDistribution(N,suzuki_sig00,suzuki_rho1,suzuki_mu12,suzuki_sigmu12);

  //Saving Suzuki Samples to a file
  OutFile.open("SuzukiSamples_figure6_7.csv", std::ios::out | std::ios::trunc);
  OutFile << "rho0_sigmu0"<< ';' ;
  OutFile << "rho0_sigmu1"<< ';' ;
  OutFile << "rho0_sigmu2"<< ';' ;
  OutFile << "rho1_sigmu0"<< ';' ;
  OutFile << "rho1_sigmu1"<< ';' ;
  OutFile << "rho1_sigmu2"<< "\n" ;

  for (int i = 0; i < N; i++) {
    OutLine.push_back(std::to_string(aux0[i]));
    OutLine.push_back(std::to_string(aux1[i]));
    OutLine.push_back(std::to_string(aux2[i]));
    OutLine.push_back(std::to_string(aux3[i]));
    OutLine.push_back(std::to_string(aux4[i]));
    OutLine.push_back(std::to_string(aux5[i]));
    copy(OutLine.begin(), OutLine.end() - 1, out_it);
    OutFile << OutLine.back() << '\n';
    OutLine.clear();
  }
  OutFile.close();

  //cleaning the vector
  aux0.clear();
  aux1.clear();
  aux2.clear();
  aux3.clear();
  aux4.clear();
  aux5.clear();

  //Deleting the channel object and finishing the program.
  delete channel;
  return 0;
}
