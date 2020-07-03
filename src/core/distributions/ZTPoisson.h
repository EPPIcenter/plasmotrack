//
// Created by Maxwell Murphy on 6/5/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRUNCATEDPOISSON_H
#define TRANSMISSION_NETWORKS_APP_TRUNCATEDPOISSON_H

#include <boost/math/special_functions/factorials.hpp>

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include "core/datatypes/Matrix.h"

#include "core/parameters/Parameter.h"

template <int MAX_COUNT>
class ZTPoisson : public Computation<ProbabilityVector<MAX_COUNT + 1>>,
                  public Observable<ZTPoisson<MAX_COUNT>>,
                  public Cacheable<ZTPoisson<MAX_COUNT>>,
                  public Checkpointable<ZTPoisson<MAX_COUNT>, ProbabilityVector<MAX_COUNT + 1>> {
public:
    explicit ZTPoisson(Parameter<double> &mean) noexcept;

    ProbabilityVector<MAX_COUNT + 1> value() noexcept;

private:
    friend class Checkpointable<ZTPoisson<MAX_COUNT>, ProbabilityVector<MAX_COUNT + 1>>;
    friend class Cacheable<ZTPoisson<MAX_COUNT>>;

    Parameter<double> &mean_;
};

template<int MAX_COUNT>
ZTPoisson<MAX_COUNT>::ZTPoisson(Parameter<double> &mean) noexcept : mean_(mean){
    this->value_(0) = 0;
    mean_.registerCacheableCheckpointTarget(this);
    mean_.add_post_change_listener([=]() { this->setDirty(); });
    this->setDirty();
    this->value();
}

template<int MAX_COUNT>
ProbabilityVector<MAX_COUNT + 1> ZTPoisson<MAX_COUNT>::value() noexcept {
    if (this->isDirty()) {
        double denominator = 0.0;
        for (int j = 1; j < MAX_COUNT + 1; ++j) {
            this->value_(j) = std::pow(mean_.value(), j) * std::exp(-mean_.value()) / boost::math::factorial<double>(j);
            denominator += this->value_(j);
        }
        this->value_ = this->value_ / denominator;
        this->setClean();
    }
    return this->value_;
}



#endif //TRANSMISSION_NETWORKS_APP_TRUNCATEDPOISSON_H
