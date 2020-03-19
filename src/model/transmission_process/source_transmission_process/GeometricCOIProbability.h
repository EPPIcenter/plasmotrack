//
// Created by Maxwell Murphy on 2/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GEOMETRICCOIPROBABILITY_H
#define TRANSMISSION_NETWORKS_APP_GEOMETRICCOIPROBABILITY_H

#include <Eigen/Core>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/computation/Computation.h"

#include "core/datatypes/Matrix.h"

#include "core/parameters/Parameter.h"


template<int MAX_COI>
class GeometricCOIProbability : public Computation<ProbabilityVector<MAX_COI + 1>>,
                                public Observable<GeometricCOIProbability<MAX_COI>>,
                                public Cacheable<GeometricCOIProbability<MAX_COI>>,
                                public Checkpointable<GeometricCOIProbability<MAX_COI>, ProbabilityVector<
                                        MAX_COI + 1>> {

public:
    explicit GeometricCOIProbability(Parameter<double> &prob);

    ProbabilityVector<MAX_COI + 1> value() noexcept;

private:
    friend class Checkpointable<GeometricCOIProbability<MAX_COI>, ProbabilityVector<MAX_COI + 1>>;

    friend class Cacheable<GeometricCOIProbability<MAX_COI>>;

    Parameter<double> &prob_;
};

template<int MAX_COI>
GeometricCOIProbability<MAX_COI>::GeometricCOIProbability(Parameter<double> &prob) : prob_(prob) {
    this->value_(0) = 0;
    prob_.registerCacheableCheckpointTarget(*this);
    prob_.add_post_change_listener([&]() { this->setDirty(); });
}

template<int MAX_COI>
ProbabilityVector<MAX_COI + 1> GeometricCOIProbability<MAX_COI>::value() noexcept {
    if (this->isDirty()) {
        for (int j = 1; j < MAX_COI + 1; ++j) {
            this->value_(j) = pow(1 - prob_.value(), j) * (prob_.value()); // geometric distribution pmf(j)
        }
        this->value_ = this->value_ / this->value_.sum();
        this->setClean();
    }

    return this->value_;
};

#endif //TRANSMISSION_NETWORKS_APP_GEOMETRICCOIPROBABILITY_H
