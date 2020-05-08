//
// Created by Maxwell Murphy on 4/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSFORMER_H
#define TRANSMISSION_NETWORKS_APP_TRANSFORMER_H


#include <cmath>
#include <functional>

#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/parameters/Parameter.h"

template<typename Output, typename Input, typename Func>
struct Transformer : public Computation<Output>,
                     public Observable<Transformer<Output, Input, Func>>,
                     public Cacheable<Transformer<Output, Input, Func>>,
                     public Checkpointable<Transformer<Output, Input, Func>, Output> {

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


#endif //TRANSMISSION_NETWORKS_APP_TRANSFORMER_H
