//
// Created by Maxwell Murphy on 4/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOGTRANSFORMER_H
#define TRANSMISSION_NETWORKS_APP_LOGTRANSFORMER_H

#include <cmath>
#include <functional>

#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/parameters/Parameter.h"

template<typename Input>
struct LogTransformer : PartialLikelihood {

    explicit LogTransformer(Input &target) : target_(target) {
        target.add_set_dirty_listener([=, this]() { this->setDirty(); });
        target.registerCacheableCheckpointTarget(this);
    }

    double value() override {
        if (this->isDirty()) {
            this->value_ = log(target_.value());
            this->setClean();
        };
        return this->value_;
    }

    Input& target_;
};

#endif //TRANSMISSION_NETWORKS_APP_LOGTRANSFORMER_H
