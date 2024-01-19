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

#include <memory>
#include <utility>
#include <vector>


namespace transmission_nets::core::computation {

    template<typename InfectionEventImpl>
    class ObservationTimeDerivedOrdering : public Computation<std::vector<std::shared_ptr<InfectionEventImpl>>>,
                                           public abstract::Observable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<std::shared_ptr<InfectionEventImpl>>> {
        using MovedCallback = std::function<void(std::shared_ptr<InfectionEventImpl> element)>;
        using ChangedCallback = std::function<void(std::shared_ptr<InfectionEventImpl> element)>;
        CREATE_KEYED_EVENT(moved_left, std::shared_ptr<InfectionEventImpl>, MovedCallback);
        CREATE_KEYED_EVENT(moved_right, std::shared_ptr<InfectionEventImpl>, MovedCallback);
        CREATE_EVENT(element_changed, ChangedCallback);

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
        void elementChanged(std::shared_ptr<InfectionEventImpl> ref);
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
//        register_element_changed_listener_key(ref);

        ref->infectionDuration()->add_post_change_listener([=, this]() { this->infectionDurationChanged(ref); });
        ref->infectionDuration()->registerCacheableCheckpointTarget(this);

        ref->add_post_change_listener([=, this]() { this->elementChanged(ref); });
        ref->registerCacheableCheckpointTarget(this);
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElements(const std::vector<std::shared_ptr<InfectionEventImpl>>& refs) noexcept {
        for (auto ref : refs) {
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

        // Find the index of the reference in the vector.
        size_t refIdx = 0;
        for (; this->value_[refIdx] != ref; ++refIdx)
            ;

        // Find the infection time of the infection just after the reference.
        // If the reference is the last element, then the infection time is a high number.
        double rightTime = std::numeric_limits<double>::max();
        if (refIdx < this->value_.size() - 1) {
            rightTime = this->value_[refIdx + 1]->infectionTime();
        }

        // Find the infection time of the infection just before the reference.
        // If the reference is the first element, then the infection time is a low number.
        double leftTime = std::numeric_limits<double>::min();
        if (refIdx > 0) {
            leftTime = (this->value_[refIdx - 1])->infectionTime();
        }

        // If the reference infection time is less than the left time or greater than the right time, then the reference is
        // out of order.
        bool unsorted = refInfectionTime < leftTime or refInfectionTime > rightTime;

        if (unsorted) {
            this->setDirty();
            // if the reference infection time is less than the left time, then the reference is moved to the left.
            if (refInfectionTime < leftTime) {
                while (unsorted) {
                    // Make sure the reference is not the first element.
                    if (refIdx > 0) {
                        // Move the reference to the left.
                        auto& curr = this->value_[refIdx];
                        auto& prev = this->value_[refIdx - 1];
                        std::swap(curr, prev);

                        // Notify the listeners that the reference has moved to the left.
                        keyed_notify_moved_left(this->value_[refIdx], this->value_[refIdx - 1]);
                        keyed_notify_moved_right(this->value_[refIdx - 1], this->value_[refIdx]);

                        // Update the reference index.
                        --refIdx;

                        // check if the reference is sorted.
                        if (refIdx == 0 or this->value_[refIdx - 1]->infectionTime() < refInfectionTime) {
                            unsorted = false;
                        }
                    } else {
                        unsorted = false;
                    }
                }
            } else {
                while (unsorted) {
                    // Make sure the reference is not the last element.
                    if (refIdx < this->value_.size() - 1) {
                        // Move the reference to the right.
                        auto& curr = this->value_[refIdx];
                        auto& next = this->value_[refIdx + 1];
                        std::swap(curr, next);

                        // Notify the listeners that the reference has moved to the right.
                        keyed_notify_moved_left(this->value_[refIdx + 1], this->value_[refIdx]);
                        keyed_notify_moved_right(this->value_[refIdx], this->value_[refIdx + 1]);

                        // Update the reference index.
                        ++refIdx;

                        // check if the reference is sorted.
                        if (refIdx == this->value_.size() - 1 or this->value_[refIdx + 1]->infectionTime() > refInfectionTime) {
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
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::elementChanged(std::shared_ptr<InfectionEventImpl> ref) {
        this->notify_element_changed(ref);
        this->setDirty();
    }

    template<typename InfectionEventImpl>
    std::vector<std::shared_ptr<InfectionEventImpl>> ObservationTimeDerivedOrdering<InfectionEventImpl>::value() {
        this->setClean();
        return this->value_;
    }

}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
