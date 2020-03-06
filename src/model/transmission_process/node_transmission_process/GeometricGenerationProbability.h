//
// Created by Maxwell Murphy on 2/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GEOMETRICGENERATIONPROBABILITY_H
#define TRANSMISSION_NETWORKS_APP_GEOMETRICGENERATIONPROBABILITY_H

#include <Eigen/Core>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/Matrix.h"

#include "core/parameters/Parameter.h"


template<int MAX_TRANSMISSIONS>
class GeometricGenerationProbability : public Computation<ProbabilityVector<MAX_TRANSMISSIONS + 1>>,
                                       public Observable<GeometricGenerationProbability<MAX_TRANSMISSIONS>>,
                                       public Cacheable<GeometricGenerationProbability<MAX_TRANSMISSIONS>>,
                                       public Checkpointable<GeometricGenerationProbability<MAX_TRANSMISSIONS>, ProbabilityVector<MAX_TRANSMISSIONS + 1>> {

public:
    explicit GeometricGenerationProbability(Parameter<double> &prob) : prob_(prob) {
        this->value_(0) = 0;
        prob_.registerCheckpointTarget(*this);
        prob_.add_post_change_listener([&]() { this->setDirty(); });
    };

    ProbabilityVector<MAX_TRANSMISSIONS + 1> value() noexcept {
        if(this->isDirty()) {
            double denominator = 0.0;
            for (int j = 1; j < MAX_TRANSMISSIONS + 1; ++j) {
                this->value_(j) = pow(1 - prob_.value(), j) * (prob_.value()); // geometric distribution pmf(j)
                denominator += this->value_(j);
            }
            this->value_ = this->value_ / denominator;
        }
        return this->value_;
    };

private:
    friend class Checkpointable<GeometricGenerationProbability<MAX_TRANSMISSIONS>, ProbabilityVector<MAX_TRANSMISSIONS + 1>>;
    friend class Cacheable<GeometricGenerationProbability<MAX_TRANSMISSIONS>>;

    Parameter<double> &prob_;
};


#endif //TRANSMISSION_NETWORKS_APP_GEOMETRICGENERATIONPROBABILITY_H
