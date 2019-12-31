//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSMISSIONMECHANISM_H
#define TRANSMISSION_NETWORKS_APP_TRANSMISSIONMECHANISM_H

#include "LikelihoodFunction.h"
#include "TransmissionMechanismParameters.h"

class TransmissionMechanism : public LikelihoodFunction<TransmissionMechanismParameters> {
public:
    LogLik evaluate(TransmissionMechanismParameters) override;
};


#endif //TRANSMISSION_NETWORKS_APP_TRANSMISSIONMECHANISM_H
