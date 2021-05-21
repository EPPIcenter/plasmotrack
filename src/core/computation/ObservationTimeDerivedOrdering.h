//
// Created by Maxwell Murphy on 3/22/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include <vector>


namespace transmission_nets::core::computation {

    template <typename InfectionEventImpl>
    class ObservationTimeDerivedOrdering : public Computation<std::vector<InfectionEventImpl*>>,
                                           public abstract::Observable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<InfectionEventImpl*>>
    {
        using MovedCallback = std::function<void(InfectionEventImpl* element)>;
        CREATE_KEYED_EVENT(moved_left, InfectionEventImpl*, MovedCallback);
        CREATE_KEYED_EVENT(moved_right, InfectionEventImpl*, MovedCallback);

    public:
        explicit ObservationTimeDerivedOrdering() noexcept;
        explicit ObservationTimeDerivedOrdering(std::vector<InfectionEventImpl*> refs) noexcept;

        void addElement(InfectionEventImpl* ref) noexcept;
        void addElements(std::vector<InfectionEventImpl*>& refs) noexcept;
        std::vector<InfectionEventImpl *> value() override;

    protected:
        friend class abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>;
        friend class abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<InfectionEventImpl*>>;

    private:

        void infectionDurationChanged(InfectionEventImpl* ref);
    };

    template<typename InfectionEventImpl>
    ObservationTimeDerivedOrdering<InfectionEventImpl>::ObservationTimeDerivedOrdering() noexcept : Computation<std::vector<InfectionEventImpl*>>() {
    }

    template<typename InfectionEventImpl>
    ObservationTimeDerivedOrdering<InfectionEventImpl>::ObservationTimeDerivedOrdering(std::vector<InfectionEventImpl *> refs) noexcept {
        addElements(refs);
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElement(InfectionEventImpl *ref) noexcept {
       this->value_.push_back(ref);
       register_moved_left_listener_key(ref);
       register_moved_right_listener_key(ref);

       ref->infectionDuration().add_post_change_listener([=, this]() { this->infectionDurationChanged(ref); });
       ref->infectionDuration().registerCacheableCheckpointTarget(this);
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElements(std::vector<InfectionEventImpl *>& refs) noexcept {
        for (auto& ref: refs) {
            addElement(ref);
        }

        std::sort(this->value_.begin(),
                  this->value_.end(),
                  [](InfectionEventImpl* a, InfectionEventImpl* b) { return a->infectionTime() < b->infectionTime(); });

        this->setDirty();
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::infectionDurationChanged(InfectionEventImpl *ref) {
        double refInfectionTime = ref->infectionTime();
        bool sourceAfterDest = false;

        auto destIt = this->value_.begin();
        auto sourceIt = destIt;

        // find where we're going to place the value -- maybe before or after the source location.
        InfectionEventImpl* destEl = *destIt;
        InfectionEventImpl* sourceEl = destEl;

        auto destInfectionTime = destEl->infectionTime();
        while (destInfectionTime < refInfectionTime) {
            destEl = *(++destIt);
            destInfectionTime = destEl->infectionTime();
            assert(destIt != this->value().end());
            if (sourceEl != ref) {
                sourceEl = destEl;
                ++sourceIt;
            }
        }

        // if the destination is before the source, we need to continue to find the source
        while(sourceEl != ref) {
           sourceEl = *(++sourceIt);
           sourceAfterDest = true;
        }


       auto betweenIt = sourceIt;
       InfectionEventImpl* betweenEl;
       // if source is not same as dest then we update and notify
       if(sourceIt != destIt) {
           // source is after dest so decrement the between iterator,
           if (sourceAfterDest) {
               do {
                   betweenEl = *(--betweenIt);
                   // betweenIt will be moved to the right of sourceIt
                   keyed_notify_moved_right(sourceEl, betweenEl);

                   // sourceIt will be moved to the left of betweenIt
                   keyed_notify_moved_left(betweenEl, sourceEl);
               }
               while (betweenIt != destIt);
           } else {
               do {
                   betweenEl = *(++betweenIt);
                   // betweenIt will be moved to the left of sourceIt
                   keyed_notify_moved_left(sourceEl, betweenEl);

                   // sourceIt will be moved to the right of betweenIt
                   keyed_notify_moved_right(betweenEl, sourceEl);
               }
               while (betweenIt != destIt);
           }
           // everything has been notified, now we erase and insert
           this->value_.erase(sourceIt);
           this->value_.insert(destIt, ref);
       }

    }
    template<typename InfectionEventImpl>
    std::vector<InfectionEventImpl *> ObservationTimeDerivedOrdering<InfectionEventImpl>::value() {
        return this->value_;
    }

}


#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
