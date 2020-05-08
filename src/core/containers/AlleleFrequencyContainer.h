//
// Created by Maxwell Murphy on 3/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H

#include <boost/container/flat_map.hpp>
#include <ostream>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"

#include "core/containers/Locus.h"

#include "core/parameters/Parameter.h"

// TODO: This container could be generalized, very similar to infection container.

template<typename AlleleFrequencyImpl, typename LocusImpl = Locus>
class AlleleFrequencyContainer : public Observable<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                 public UncacheablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                 public CheckpointablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>> {

    using AlleleFrequencyMap = boost::container::flat_map<LocusImpl *, Parameter<AlleleFrequencyImpl>>;

public:
    using LocusAlleleFrequencyAssignment = std::pair<LocusImpl *, AlleleFrequencyImpl>;

    std::vector<LocusImpl *> loci{};

    AlleleFrequencyContainer() = default;

    template<typename AlleleDataIter>
    explicit AlleleFrequencyContainer(AlleleDataIter alleleFreqs);

    AlleleFrequencyMap& alleleFrequencies() {
        return alleleFrequencies_;
    };

    Parameter<AlleleFrequencyImpl>& alleleFrequencies(LocusImpl &locus) {
        return alleleFrequencies_.at(&locus);
    };

    Parameter<AlleleFrequencyImpl>& alleleFrequencies(LocusImpl *locus) {
        return alleleFrequencies_.at(locus);
    };


    void addLocus(LocusImpl &locus);

    void addLocus(LocusImpl *locus);

    friend std::ostream &operator<<(std::ostream &os, const AlleleFrequencyContainer &container) {
        for (const auto& [locus, alleleFreqs] : container.alleleFrequencies_ ) {
            os << locus->label << ": " << alleleFreqs.value() << std::endl;
        }
        return os;
    };


private:
    AlleleFrequencyMap alleleFrequencies_{};

};

template<typename AlleleFrequencyImpl, typename LocusImpl>
template<typename AlleleDataIter>
AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::AlleleFrequencyContainer(AlleleDataIter alleleFreqs) {
    for (auto &[locus, alleleFrequencies] : alleleFreqs) {
        assert(locus->totalAlleles() == alleleFrequencies.totalElements());
        loci.push_back(locus);
        alleleFrequencies_.emplace(locus, alleleFrequencies);
        // pass through notifications
//        alleleFrequencies_.at(locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
        alleleFrequencies_.at(locus).add_post_change_listener([=]() { this->notify_post_change(); });
        alleleFrequencies_.at(locus).add_save_state_listener([=]() { this->notify_save_state(); });
        alleleFrequencies_.at(locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
        alleleFrequencies_.at(locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
    }
}

template<typename AlleleFrequencyImpl, typename LocusImpl>
void AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::addLocus(LocusImpl &locus) {
    loci.push_back(&locus);
    alleleFrequencies_.emplace(&locus, AlleleFrequencyImpl(locus.totalAlleles()));
//    alleleFrequencies_.at(&locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
    alleleFrequencies_.at(&locus).add_post_change_listener([=]() {
        this->notify_post_change();
    });
    alleleFrequencies_.at(&locus).add_save_state_listener([=]() { this->notify_save_state(); });
    alleleFrequencies_.at(&locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
    alleleFrequencies_.at(&locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
}

template<typename AlleleFrequencyImpl, typename LocusImpl>
void AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::addLocus(LocusImpl *locus) {
    loci.push_back(locus);
    alleleFrequencies_.emplace(locus, AlleleFrequencyImpl(locus->totalAlleles()));
//    alleleFrequencies_.at(locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
    alleleFrequencies_.at(locus).add_post_change_listener([=]() {
        this->notify_post_change();
    });
    alleleFrequencies_.at(locus).add_save_state_listener([=]() { this->notify_save_state(); });
    alleleFrequencies_.at(locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
    alleleFrequencies_.at(locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
}

//template<typename AlleleFrequencyImpl, typename LocusImpl>
//void AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::notify_post_change() const {
//    UncacheablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>::notify_post_change();
//}

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
