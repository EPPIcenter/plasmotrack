//
// Created by Maxwell Murphy on 3/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H

#include <boost/container/flat_map.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Matrix.h"
#include "core/datatypes/AlleleFrequenciesVector.h"
#include "core/parameters/Parameter.h"

// TODO: This container could be generalized, very similar to infection container.

template<typename AlleleFrequencyImpl, typename LocusImpl = Locus>
class AlleleFrequencyContainer : public Observable<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                 public UncacheablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                 public CheckpointablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>> {

    using AlleleFrequencyMap = boost::container::flat_map<LocusImpl *, Parameter<AlleleFrequencyImpl>>;

public:
    using LocusAlleleFrequencyAssignment = std::pair<LocusImpl *, AlleleFrequencyImpl>;

    template<typename AlleleDataIter>
    explicit AlleleFrequencyContainer(AlleleDataIter alleleFreqs);

    AlleleFrequencyMap alleleFrequencies() {
        return alleleFrequencies_;
    };

    Parameter<AlleleFrequencyImpl> alleleFrequencies(LocusImpl *locus) {
        return alleleFrequencies_.at(locus);
    };

    Parameter<AlleleFrequencyImpl> alleleFrequencies(LocusImpl *locus) const {
        return alleleFrequencies_.at(locus);
    };


private:
    AlleleFrequencyMap alleleFrequencies_{};
};

template<typename AlleleFrequencyImpl, typename LocusImpl>
template<typename AlleleDataIter>
AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::AlleleFrequencyContainer(AlleleDataIter alleleFreqs) {
    for (auto const&[locus, alleleFrequencies] : alleleFreqs) {
        assert(locus->totalAlleles() == alleleFrequencies.totalAlleles());
        alleleFrequencies_.emplace(locus, alleleFrequencies);
        // pass through notifications
        alleleFrequencies_.at(locus).add_pre_change_listener([&]() { this->notify_pre_change(); });
        alleleFrequencies_.at(locus).add_post_change_listener([&]() { this->notify_post_change(); });
        alleleFrequencies_.at(locus).add_save_state_listener([&]() { this->notify_save_state(); });
        alleleFrequencies_.at(locus).add_accept_state_listener([&]() { this->notify_accept_state(); });
        alleleFrequencies_.at(locus).add_restore_state_listener([&]() { this->notify_restore_state(); });
    }
}

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
