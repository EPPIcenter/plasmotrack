//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
#define TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H

#include <boost/container/flat_set.hpp>

#include "Computation.h"
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

    void customSetClean() noexcept;

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
    dirty_targets_.insert(&target);

    target.add_set_dirty_listener([&]() {
        dirty_targets_.insert(&target);
        this->value_ -= target.peek();
        this->setDirty();
    });

    target.registerCheckpointTarget(*this);
}

template<typename Input, typename Output>
void Accumulator<Input, Output>::customSetClean() noexcept {
    this->dirty_targets_.clear();
    this->setClean();
}

template<typename Input, typename Output>
Output Accumulator<Input, Output>::value() noexcept {
    for (auto el : dirty_targets_) {
        this->value_ += el->value();
    }
    customSetClean();
    return this->value_;
}

template<typename Input, typename Output>
int Accumulator<Input, Output>::getNumTargets() const {
    return targets_.size();
}

#endif //TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
