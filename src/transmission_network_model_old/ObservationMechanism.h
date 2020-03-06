//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H

#include "LikelihoodFunction.h"
#include "ObservationMechanismParameters.h"

class ObservationMechanism : public LikelihoodFunction<ObservationMechanismParameters> {
public:
    LogLik evaluate(ObservationMechanismParameters params) override;
};


#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H
