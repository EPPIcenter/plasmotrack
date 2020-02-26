//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include "core/containers/Infection.h"

template<int MAX_COI>
class MultinomialSourceTransmissionProcess : public Observable<MultinomialSourceTransmissionProcess<MAX_COI>>,
                                             public Cacheable<MultinomialSourceTransmissionProcess<MAX_COI>> {

    //  TODO: COI Probability -- need to implement a distribution over COI instead of just flat I think.
public:
    template<typename GeneticsImpl>
    double calculateLogLikelihood(Infection<GeneticsImpl>& founder) {
        return 0;
    }


private:

};

#endif //TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
