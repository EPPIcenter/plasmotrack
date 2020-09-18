//
// Created by Maxwell Murphy on 4/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSFORMER_H
#define TRANSMISSION_NETWORKS_APP_TRANSFORMER_H

#include "core/computation/Computation.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

namespace transmission_nets::core::computation {

    template<typename Output, typename Input, typename Func>
    struct Transformer : public Computation<Output>,
                         public abstract::Observable<Transformer<Output, Input, Func>>,
                         public abstract::Cacheable<Transformer<Output, Input, Func>>,
                         public abstract::Checkpointable<Transformer<Output, Input, Func>, Output> {

        Transformer(Input &target, Func tFunc) : target_(target), tFunc_(tFunc) {}

        Output value() override {
            if (this->isDirty()) {
                this->value_ = tfunc(target_.value());
                this->setClean();
            };
            return this->value_;
        }

        Input& target_;
        Func& tFunc_;
    };

}



#endif //TRANSMISSION_NETWORKS_APP_TRANSFORMER_H
