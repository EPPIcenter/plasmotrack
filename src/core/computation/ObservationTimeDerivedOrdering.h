//
// Created by Maxwell Murphy on 3/22/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include <fmt/core.h>

#include <vector>
#include <utility>
#include <memory>


namespace transmission_nets::core::computation {

    template <typename InfectionEventImpl>
    class ObservationTimeDerivedOrdering : public Computation<std::vector<std::shared_ptr<InfectionEventImpl>>>,
                                           public abstract::Observable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<std::shared_ptr<InfectionEventImpl>>>
    {
        using MovedCallback = std::function<void(std::shared_ptr<InfectionEventImpl> element)>;
        CREATE_KEYED_EVENT(moved_left, std::shared_ptr<InfectionEventImpl>, MovedCallback);
        CREATE_KEYED_EVENT(moved_right, std::shared_ptr<InfectionEventImpl>, MovedCallback);

    public:
        explicit ObservationTimeDerivedOrdering() noexcept;
        explicit ObservationTimeDerivedOrdering(const std::vector<std::shared_ptr<InfectionEventImpl>>& refs) noexcept;

        void addElements(const std::vector<std::shared_ptr<InfectionEventImpl>>& refs) noexcept;
        std::vector<std::shared_ptr<InfectionEventImpl>> value() override;

    protected:
        friend class abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>;
        friend class abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<std::shared_ptr<InfectionEventImpl>>>;

    private:
        void infectionDurationChanged(std::shared_ptr<InfectionEventImpl> ref);
        void addElement(std::shared_ptr<InfectionEventImpl> ref) noexcept;
    };

    template<typename InfectionEventImpl>
    ObservationTimeDerivedOrdering<InfectionEventImpl>::ObservationTimeDerivedOrdering() noexcept : Computation<std::vector<std::shared_ptr<InfectionEventImpl>>>() {
    }

    template<typename InfectionEventImpl>
    ObservationTimeDerivedOrdering<InfectionEventImpl>::ObservationTimeDerivedOrdering(const std::vector<std::shared_ptr<InfectionEventImpl>>& refs) noexcept {
        addElements(refs);
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElement(std::shared_ptr<InfectionEventImpl> ref) noexcept {
       this->value_.push_back(ref);
       register_moved_left_listener_key(ref);
       register_moved_right_listener_key(ref);

       ref->infectionDuration()->add_post_change_listener([=, this]() { this->infectionDurationChanged(ref); });
       ref->infectionDuration()->registerCacheableCheckpointTarget(this);
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElements(const std::vector<std::shared_ptr<InfectionEventImpl>>& refs) noexcept {
        for (auto ref: refs) {
            addElement(std::move(ref));
        }

        std::sort(this->value_.begin(),
                  this->value_.end(),
                  [](std::shared_ptr<InfectionEventImpl> a, std::shared_ptr<InfectionEventImpl> b) { return a->infectionTime() < b->infectionTime(); });

        this->setDirty();
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::infectionDurationChanged(std::shared_ptr<InfectionEventImpl> ref) {
        double refInfectionTime = ref->infectionTime();

        size_t refIdx = 0;
        for(; this->value_[refIdx] != ref; ++refIdx);

        double rightTime = std::numeric_limits<double>::max();
        if (refIdx < this->value_.size() - 1) {
            rightTime = this->value_[refIdx + 1]->infectionTime();
        }

        double leftTime = std::numeric_limits<double>::min();
        if (refIdx > 0) {
            leftTime = (this->value_[refIdx - 1])->infectionTime();
        }

        bool unsorted = refInfectionTime < leftTime or refInfectionTime > rightTime;

        if (unsorted) {
            this->setDirty();
            if (refInfectionTime < leftTime){
                while(unsorted) {
                    if (refIdx > 0) {
                        auto& curr = this->value_[refIdx];
                        auto& prev = this->value_[refIdx - 1];
                        std::swap(curr, prev);
                        keyed_notify_moved_left(this->value_[refIdx], this->value_[refIdx - 1]);
                        keyed_notify_moved_right(this->value_[refIdx - 1], this->value_[refIdx]);
                        --refIdx;
                        if(refIdx == 0 or this->value_[refIdx - 1]->infectionTime() < refInfectionTime) {
                            unsorted = false;
                        }
                    } else {
                        unsorted = false;
                    }
                }
            } else {
                while(unsorted) {
                    if (refIdx < this->value_.size() - 1) {
                        auto& curr = this->value_[refIdx];
                        auto& next = this->value_[refIdx + 1];
                        std::swap(curr, next);
                        keyed_notify_moved_left(this->value_[refIdx + 1], this->value_[refIdx]);
                        keyed_notify_moved_right(this->value_[refIdx], this->value_[refIdx + 1]);
                        ++refIdx;
                        if(refIdx == this->value_.size() - 1 or this->value_[refIdx + 1]->infectionTime() > refInfectionTime) {
                            unsorted = false;
                        }
                    } else {
                        unsorted = false;
                    }
               }
            }
        }
    }

    template<typename InfectionEventImpl>
    std::vector<std::shared_ptr<InfectionEventImpl>> ObservationTimeDerivedOrdering<InfectionEventImpl>::value() {
        this->setClean();
        return this->value_;
    }

}


#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
