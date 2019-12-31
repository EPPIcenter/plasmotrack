//
// Created by Maxwell Murphy on 12/10/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H

#include "PartialLikelihood.h"
#include "../Parameters/SimpleParameter.h"

class AdderLikelihood : public PartialLikelihood {
private:
    IntegerParameter* a_;
    IntegerParameter* b_;
public:
    AdderLikelihood(IntegerParameter *a, IntegerParameter *b);

    float value() override;

    float peek() override;
};


#endif //TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H
