/** Implementation of the Process Elementary Properties class */

#ifndef PROCESSELEMETARYPROPERTIES_HPP
#define PROCESSELEMETARYPROPERTIES_HPP

#include <array>

template<short N, typename fpt>

class ProcessElementaryProperties{
    public:
    fpt CalcMeanValue();
    fpt CalcMeanPower(std::array<fpt,N> &pathGains);
  };

//document me
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanValue(){return 0;}


//document me
template<short N, typename fpt>
fpt ProcessElementaryProperties<N,fpt>::CalcMeanPower(std::array<fpt,N> &pathGains){
    fpt MeanPower=0;
    for(int i=0;i<N;i++){ 
      MeanPower+=pow(pathGains[i],2)/2;
    }
    return MeanPower;
}


#endif
