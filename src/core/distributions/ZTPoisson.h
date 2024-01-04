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

namespace transmission_nets::core::distributions {
    template<int MAX_COUNT>
    class ZTPoisson : public computation::Computation<std::array<long double, MAX_COUNT + 1>>,
                      public abstract::Observable<ZTPoisson<MAX_COUNT>>,
                      public abstract::Cacheable<ZTPoisson<MAX_COUNT>>,
                      public abstract::Checkpointable<ZTPoisson<MAX_COUNT>, std::array<long double, MAX_COUNT + 1>> {
        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;

    public:
        explicit ZTPoisson(p_ParameterDouble mean) noexcept;

        std::array<long double, MAX_COUNT + 1> value() noexcept;

    private:
        friend class abstract::Checkpointable<ZTPoisson<MAX_COUNT>, std::array<long double, MAX_COUNT + 1>>;
        friend class abstract::Cacheable<ZTPoisson<MAX_COUNT>>;

        p_ParameterDouble mean_;
    };

    template<int MAX_COUNT>
    ZTPoisson<MAX_COUNT>::ZTPoisson(p_ParameterDouble mean) noexcept : mean_(std::move(mean)) {
        this->value_[0] = -std::numeric_limits<long double>::infinity();
        mean_->registerCacheableCheckpointTarget(this);
        mean_->add_post_change_listener([=, this]() {
            this->setDirty();
        });
        this->setDirty();
        this->value();
    }

    template<int MAX_COUNT>
    std::array<long double, MAX_COUNT + 1> ZTPoisson<MAX_COUNT>::value() noexcept {
        if (this->isDirty()) {
            double lambda = mean_->value();
            long double denominator = 0.0;
            for (int j = 1; j < MAX_COUNT + 1; ++j) {
                this->value_[j] = j * std::log(lambda) - std::log(std::exp(lambda) - 1) - std::log(boost::math::factorial<double>(j));
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

}// namespace transmission_nets::core::distributions


#endif//TRANSMISSION_NETWORKS_APP_TRUNCATEDPOISSON_H
