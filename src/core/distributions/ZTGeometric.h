//
// Created by Maxwell Murphy on 2/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZTGEOMETRIC_H
#define TRANSMISSION_NETWORKS_APP_ZTGEOMETRIC_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include "core/datatypes/Matrix.h"

#include "core/parameters/Parameter.h"

#include <memory>

namespace transmission_nets::core::distributions {

    using namespace transmission_nets::core::parameters;

    template<int MAX_COUNT>
    class ZTGeometric : public computation::Computation<datatypes::ProbabilityVector<MAX_COUNT + 1>>,
                        public abstract::Observable<ZTGeometric<MAX_COUNT>>,
                        public abstract::Cacheable<ZTGeometric<MAX_COUNT>>,
                        public abstract::Checkpointable<ZTGeometric<MAX_COUNT>, datatypes::ProbabilityVector<
                                                                                        MAX_COUNT + 1>> {

    public:
        explicit ZTGeometric(std::shared_ptr<Parameter<double>> prob) noexcept;

        datatypes::ProbabilityVector<MAX_COUNT + 1> value() noexcept override;

    private:
        friend class abstract::Checkpointable<ZTGeometric, datatypes::ProbabilityVector<
                                                                              MAX_COUNT + 1>>;

        friend class abstract::Cacheable<ZTGeometric>;

        std::shared_ptr<Parameter<double>> prob_;
    };

    template<int MAX_COUNT>
    ZTGeometric<MAX_COUNT>::ZTGeometric(std::shared_ptr<Parameter<double>> prob) noexcept : prob_(std::move(prob)) {
        this->value_(0) = 0;
        prob_->registerCacheableCheckpointTarget(this);
        prob_->add_post_change_listener([=, this]() { this->setDirty(); });
        this->setDirty();
        this->ZTGeometric::value();
    }


    template<int MAX_COUNT>
    datatypes::ProbabilityVector<MAX_COUNT + 1> ZTGeometric<MAX_COUNT>::value() noexcept {
        if (this->isDirty()) {
            double denominator = 0.0;
            for (int j = 1; j < MAX_COUNT + 1; ++j) {
                // geometric distribution pmf(j)
                this->value_(j) = pow(1 - prob_->value(), j) * (prob_->value());
                denominator += this->value_(j);
            }
            this->value_ = this->value_ / denominator;
            this->setClean();
        }
        return this->value_;
    }

}// namespace transmission_nets::core::distributions


#endif//TRANSMISSION_NETWORKS_APP_ZTGEOMETRIC_H
