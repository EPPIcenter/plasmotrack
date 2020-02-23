//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H

#include "core/computation/Computation.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

using Likelihood = double;

class PartialLikelihood : public Computation<Likelihood>,
                          public Observable<PartialLikelihood>,
                          public Cacheable<PartialLikelihood>,
                          public Checkpointable<PartialLikelihood, Likelihood> {
public:

    explicit PartialLikelihood(std::string id);

    Likelihood value() override = 0;

protected:
    friend class Cacheable<PartialLikelihood>;
    friend class Checkpointable<PartialLikelihood, Likelihood>;

    const std::string id_;
};


#endif //TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H
