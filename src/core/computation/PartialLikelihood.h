//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H

#include "core/computation/Computation.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

namespace transmission_nets::core::computation {

    using Likelihood = double;

    class PartialLikelihood : public Computation<Likelihood>,
                              public abstract::Observable<PartialLikelihood>,
                              public abstract::Cacheable<PartialLikelihood>,
                              public abstract::Checkpointable<PartialLikelihood, Likelihood> {
    public:

//    explicit PartialLikelihood() = default;

        virtual Likelihood value() override = 0;
        virtual std::string identifier() = 0;

    protected:
        friend class abstract::Cacheable<PartialLikelihood>;
        friend class abstract::Checkpointable<PartialLikelihood, Likelihood>;
    };
}



#endif //TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H
