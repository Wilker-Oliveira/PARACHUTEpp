#include "SoSModel.hpp"

template<short N, typename fpt>

class MEAModel: public SoSModel<N,fpt>{

public:

    MEAModel(float sig, float fmax){
        DefineModel(sig, fmax);
        this->genPhases();
    }

    MEAModel(float sig, float fc, float kc){
        DefineModel(sig, fc, kc);
        this->genPhases();
    }

    double rootFunc(int n, float fc, float TOL = 0.05) {

        auto f = [fc](auto dfreq){return (n/N) - erf(dfreq*sqrt(log(2))/fc)};

        // preciso definir o intervalo e tolerância
        float a,b;
        float x;

        while (diff > TOL){
            diff = fabs(a-b)/2;
            x = (a+b)/2;

            if(f(x) == 0){
                return x;
                break;
            }
            if(f(a)*f(x) < 0){
                b = x;
            }
            else a = x;
            }

        return x;
    }
};

