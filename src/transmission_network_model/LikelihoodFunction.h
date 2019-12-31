//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_LIKELIHOODFUNCTION_H
#define TRANSMISSION_NETWORKS_APP_LIKELIHOODFUNCTION_H

#include "LikelihoodFunctionParameters.h"

typedef float LogLik;

//template <class LikelihoodFunctionParameters>
class LikelihoodFunction {
public:
    virtual LogLik evaluate(LikelihoodFunctionParameters params) = 0;
};


#endif //TRANSMISSION_NETWORKS_APP_LIKELIHOODFUNCTION_H
