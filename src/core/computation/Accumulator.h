//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
#define TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H


#include "Computation.h"
#include "PartialLikelihood.h"

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Dependency.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <memory>


namespace transmission_nets::core::computation {

    template<typename Input, typename Output>
    class Accumulator : public Computation<Output>,
                        public abstract::Observable<Accumulator<Input, Output>>,
                        public abstract::Cacheable<Accumulator<Input, Output>>,
                        public abstract::Checkpointable<Accumulator<Input, Output>, Output> {

    public:
        Accumulator();
        void addTarget(const std::shared_ptr<Input>& target);

        Output value() noexcept override;

        [[nodiscard]] int getNumTargets() const;

    private:
        friend class abstract::Cacheable<Accumulator>;

        friend class abstract::Checkpointable<Accumulator, Output>;

        using TargetSet = std::vector<std::shared_ptr<Input>>;

        TargetSet targets_{};
        TargetSet dirtyTargets_{};

        // Takes `target`'s contribution out of the running sum the first time it goes dirty, until value() adds its
        // new one back.
        void targetWentDirty(const std::shared_ptr<Input>& target);
    };

    template<typename Input, typename Output>
    void Accumulator<Input, Output>::targetWentDirty(const std::shared_ptr<Input>& target) {
        for (auto& el : dirtyTargets_) {
            if (el == target) {
                return;
            }
        }
        dirtyTargets_.emplace_back(target);
        this->value_ -= target->peek();
    }


    template<typename Input, typename Output>
    void Accumulator<Input, Output>::addTarget(const std::shared_ptr<Input>& target) {
        this->setDirty();
        targets_.emplace_back(target);
        dirtyTargets_.emplace_back(target);
        abstract::dependsOn(this, target, [=, this]() { targetWentDirty(target); });
    }


    template<>
    inline void Accumulator<PartialLikelihood, Likelihood>::addTarget(const std::shared_ptr<PartialLikelihood>& target) {
        this->setDirty();

        //check if target is already in targets_
        bool inserted = false;
        for (auto& el : targets_) {
            if (el == target) {
                inserted = true;
                break;
            }
        }
#ifndef DNDEBUG
        if (inserted) assert(!"Target added more than once. Check model specification.");
#endif
        targets_.emplace_back(target);
        dirtyTargets_.emplace_back(target);
        abstract::dependsOn(this, target, [=, this]() { targetWentDirty(target); });
    }


    template<typename Input, typename Output>
    Output Accumulator<Input, Output>::value() noexcept {
        for (auto el : dirtyTargets_) {
            this->value_ += el->value();
        }
        this->setClean();
        dirtyTargets_.resize(0);

        return this->value_;
    }


    template<>
    inline Likelihood Accumulator<PartialLikelihood, Likelihood>::value() noexcept {
        assert(this->value_ < std::numeric_limits<Likelihood>::infinity());

        for (auto& el : dirtyTargets_) {
            assert(this->value_ < std::numeric_limits<Likelihood>::infinity());

            if (std::isnan(el->value())) {
                this->value_ = -std::numeric_limits<Likelihood>::infinity();
            } else {
                this->value_ += el->value();
            }
            assert(this->value_ < std::numeric_limits<Likelihood>::infinity());
        }
        this->setClean();
        dirtyTargets_.clear();

        return this->value_;
    }


    template<typename Input, typename Output>
    int Accumulator<Input, Output>::getNumTargets() const {
        return targets_.size();
    }


    template<typename Input, typename Output>
    Accumulator<Input, Output>::Accumulator() = default;

}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_ACCUMULATOR_H
