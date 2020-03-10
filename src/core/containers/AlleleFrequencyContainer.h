//
// Created by Maxwell Murphy on 3/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H

#include <boost/container/flat_map.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Matrix.h"
#include "core/parameters/Parameter.h"

// TODO: This container could be generalized, very similar to infection container.
template<typename AlleleFrequencyImpl, typename LocusImpl = Locus>
using LocusAlleleFrequencyAssignment = std::pair<LocusImpl*, AlleleFrequencyImpl>;

template<typename AlleleFrequencyImpl, typename LocusImpl = Locus>
class AlleleFrequencyContainer : public Observable<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>> {
    using AlleleFrequencyMap = boost::container::flat_map<LocusImpl*, Parameter<AlleleFrequencyImpl>>;

    using ChangeCallback = std::function<void(LocusImpl*)>;
    CREATE_EVENT(pre_change, ChangeCallback);
    CREATE_EVENT(post_change, ChangeCallback);
    using CallbackType = std::function<void()>;
    CREATE_EVENT(save_state, CallbackType);
    CREATE_EVENT(accept_state, CallbackType);
    CREATE_EVENT(restore_state, CallbackType);

public:
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0& target) {
        ListenerId_t saveStateEventId = this->add_save_state_listener([&]() { target.saveState(); });
        ListenerId_t acceptStateEventId = this->add_accept_state_listener([&]() { target.acceptState(); });
        ListenerId_t restoreStateEventId = this->add_restore_state_listener([&]() { target.restoreState(); });
        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename AlleleDataIter>
    AlleleFrequencyContainer(AlleleDataIter alleleFreqs) {
        for(auto const& [locus, alleleFrequencies] : alleleFreqs) {
            alleleFrequencies_.emplace(locus, alleleFrequencies);
            alleleFrequencies_.at(locus).add_pre_change_listener([&]() { notify_pre_change(); });
            alleleFrequencies_.at(locus).add_post_change_listener([&]() { notify_post_change(); });
            alleleFrequencies_.at(locus).add_save_state_listener([&]() { notify_save_state(); });
            alleleFrequencies_.at(locus).add_accept_state_listener([&]() { notify_accept_state(); });
            alleleFrequencies_.at(locus).add_restore_state_listener([&]() { notify_restore_state(); });
        }
    };

    AlleleFrequencyMap alleleFrequencies() {
        return alleleFrequencies_;
    };

    Parameter<AlleleFrequencyImpl> alleleFrequencies(LocusImpl* locus) {
        return alleleFrequencies_.at(locus);
    };


private:
    AlleleFrequencyMap alleleFrequencies_{};
};

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
