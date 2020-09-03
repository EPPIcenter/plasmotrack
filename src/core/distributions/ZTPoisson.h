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

#include "core/utils/numerics.h"

template <int MAX_COUNT>
class ZTPoisson : public Computation<std::array<long double, MAX_COUNT + 1>>,
                  public Observable<ZTPoisson<MAX_COUNT>>,
                  public Cacheable<ZTPoisson<MAX_COUNT>>,
                  public Checkpointable<ZTPoisson<MAX_COUNT>, std::array<long double, MAX_COUNT + 1>> {
public:
    explicit ZTPoisson(Parameter<double> &mean) noexcept;

//    ProbabilityVector<MAX_COUNT + 1> value() noexcept;
    std::array<long double, MAX_COUNT + 1> value() noexcept;

private:
    friend class Checkpointable<ZTPoisson<MAX_COUNT>,std::array<long double, MAX_COUNT + 1>>;
    friend class Cacheable<ZTPoisson<MAX_COUNT>>;

    Parameter<double> &mean_;
};

template<int MAX_COUNT>
ZTPoisson<MAX_COUNT>::ZTPoisson(Parameter<double> &mean) noexcept : mean_(mean){
    this->value_[0] = -std::numeric_limits<long double>::infinity();
    mean_.registerCacheableCheckpointTarget(this);
    mean_.add_post_change_listener([=, this]() { this->setDirty(); });
    this->setDirty();
    this->value();
}

template<int MAX_COUNT>
std::array<long double, MAX_COUNT + 1> ZTPoisson<MAX_COUNT>::value() noexcept {
    if (this->isDirty()) {

        long double denominator = 0.0;
        for (int j = 1; j < MAX_COUNT + 1; ++j) {
            this->value_[j] = j * log(mean_.value()) - mean_.value() - log(boost::math::factorial<double>(j));
            denominator += exp(this->value_[j]);
            assert(!std::isnan(this->value_[j]));
        }

        denominator = log(denominator);
        for (int i = 1; i < MAX_COUNT + 1; ++i) {
            this->value_[i] -= denominator;
            assert(!std::isnan(this->value_[i]));
        }
        this->setClean();
    }
    return this->value_;
}



#endif //TRANSMISSION_NETWORKS_APP_TRUNCATEDPOISSON_H
