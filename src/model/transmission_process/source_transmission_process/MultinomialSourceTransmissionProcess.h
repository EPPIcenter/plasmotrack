//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include <core/containers/AlleleFrequencyContainer.h>
#include "core/containers/Infection.h"

template<int MAX_COI, typename AlleleFrequencyImpl>
class MultinomialSourceTransmissionProcess : public Observable<MultinomialSourceTransmissionProcess<MAX_COI, AlleleFrequencyImpl>>,
                                             public Cacheable<MultinomialSourceTransmissionProcess<MAX_COI, AlleleFrequencyImpl>> {

    //  TODO: COI Probability -- need to implement a distribution over COI instead of just flat I think.
public:

    template<typename GeneticsImpl>
    double calculateLikelihood(Infection<GeneticsImpl>& founder) {
        return 0;
    }

    template<typename GeneticsImpl>
    double calculateLogLikelihood(Infection<GeneticsImpl>& founder) {
        return log(calculateLikelihood(founder));
    }


private:

    AlleleFrequencyContainer<AlleleFrequencyImpl> &alleleFrequencies_;
};

#endif //TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
