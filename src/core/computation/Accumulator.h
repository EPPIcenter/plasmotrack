//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
#define TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H

#include <boost/container/flat_set.hpp>
#include <cmath>

#include "Computation.h"
#include "PartialLikelihood.h"

#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Observable.h"

template<typename Input, typename Output>
class Accumulator : public Computation<Output>,
                    public Observable<Accumulator<Input, Output>>,
                    public Cacheable<Accumulator<Input, Output>>,
                    public Checkpointable<Accumulator<Input, Output>, Output> {

public:

    void addTarget(Input &target);

    void addTarget(Input *target);


    Output value() noexcept override;

    [[nodiscard]] int getNumTargets() const;

private:
    friend class Cacheable<Accumulator<Input, Output>>;

    friend class Checkpointable<Accumulator<Input, Output>, Output>;

    boost::container::flat_set<Input *> targets_{};
    boost::container::flat_set<Input *> dirty_targets_{};

};

template<typename Input, typename Output>
void Accumulator<Input, Output>::addTarget(Input &target) {
    targets_.insert(&target);
//    this->value_ += target.value();
    dirty_targets_.insert(&target);

    target.add_set_dirty_listener([&]() {
        const auto& [_, inserted] = dirty_targets_.insert(&target);
        if (inserted) {
            this->value_ -= target.peek();
            this->setDirty();
        }
    });

    target.registerCacheableCheckpointTarget(this);
    this->setDirty();
}

template<>
inline void Accumulator<PartialLikelihood, double>::addTarget(PartialLikelihood &target) {
    const auto& [_, inserted] = targets_.insert(&target);
    if(!inserted) assert(!"Target added more than once. Check model specification.");
    this->value_ += target.value();
//    dirty_targets_.insert(&target);

    target.add_set_dirty_listener([&]() {
        const auto& [_, inserted] = dirty_targets_.insert(&target);
        if (inserted) {
            this->value_ -= target.peek();
            this->setDirty();
        }
    });

    target.registerCacheableCheckpointTarget(this);
//    this->setDirty();
}

template<typename Input, typename Output>
void Accumulator<Input, Output>::addTarget(Input *target) {
    targets_.insert(target);
    this->value_ += target->value();
//    dirty_targets_.insert(target);

    target->add_set_dirty_listener([=]() {
        const auto& [_, inserted] = dirty_targets_.insert(target);
        if (inserted) {
            this->value_ -= target->peek();
            this->setDirty();
        }
    });

    target->registerCacheableCheckpointTarget(this);
//    this->setDirty();
}

template<>
inline void Accumulator<PartialLikelihood, double>::addTarget(PartialLikelihood *target) {
    const auto& [_, inserted] = targets_.insert(target);
    if(!inserted) assert(!"Target added more than once. Check model specification.");
    this->value_ += target->value();
//    dirty_targets_.insert(target);

    target->add_set_dirty_listener([=]() {
        const auto& [_, inserted] = dirty_targets_.insert(target);
        if (inserted) {
            this->setDirty();
            this->value_ -= target->peek();
        }
    });

    target->registerCacheableCheckpointTarget(this);
//    this->setDirty();
}

//template<typename Input, typename Output>
//void Accumulator<Input, Output>::customSetClean() noexcept {
//    this->dirty_targets_.clear();
//    this->setClean();
//}

template<typename Input, typename Output>
Output Accumulator<Input, Output>::value() noexcept {
    for (auto el : dirty_targets_) {
        this->value_ += el->value();
    }
    this->dirty_targets_.clear();

    this->setClean();
    return this->value_;
}


template<>
inline double Accumulator<PartialLikelihood, double>::value() noexcept {
    for (auto &el : dirty_targets_) {
        this->value_ += std::isnan(el->value()) ? std::numeric_limits<double>::lowest() : el->value();
    }
    this->dirty_targets_.clear();

    this->setClean();
    return this->value_;
}

template<typename Input, typename Output>
int Accumulator<Input, Output>::getNumTargets() const {
    return targets_.size();
}

#endif //TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
