//
// Created by Maxwell Murphy on 1/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFECTION_H
#define TRANSMISSION_NETWORKS_APP_INFECTION_H

#include <boost/container/flat_map.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Data.h"
#include "core/parameters/Parameter.h"



template<typename GeneticImpl, typename LocusImpl = Locus>
using LocusGeneticsAssignment = std::pair<LocusImpl*, GeneticImpl>;

template <typename GeneticImpl, typename LocusImpl = Locus>
class Infection: public Observable<Infection<GeneticImpl, LocusImpl>> {
    template<template <typename> typename Wrapper>
    using GenotypeMap = boost::container::flat_map<LocusImpl*, Wrapper<GeneticImpl>>; // Maps genetic readout to loci

    using ChangeCallback = std::function<void()>;
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

    template<typename LocusDataIter>
    Infection(LocusDataIter obs, LocusDataIter latent) {
        for(auto const& [locus, genetics] : obs) {
            assert(locus->totalAlleles() == genetics.totalAlleles());
            observedGenotype_.emplace(locus, genetics);
        }

        //// TODO: Generalize this notification forwarding in collections of parameters
        for(auto const& [locus, genetics] : latent) {
            latentGenotype_.emplace(locus, genetics);
            // Creating pass through of notifications
            latentGenotype_.at(locus).add_pre_change_listener([&]() { notify_pre_change(); });
            latentGenotype_.at(locus).add_post_change_listener([&]() { notify_post_change(); });
            latentGenotype_.at(locus).add_save_state_listener([&]() { notify_save_state(); });
            latentGenotype_.at(locus).add_accept_state_listener([&]() { notify_accept_state(); });
            latentGenotype_.at(locus).add_restore_state_listener([&]() { notify_restore_state(); });
        }
    };

    GenotypeMap<Data>& observedGenotype() {
        return observedGenotype_;
    };

    Data<GeneticImpl>& observedGenotype(LocusImpl* locus) {
        return observedGenotype_.at(locus);
    };

    GenotypeMap<Parameter>& latentGenotype() {
        return latentGenotype_;
    };

    Parameter<GeneticImpl>& latentGenotype(LocusImpl* locus) {
        return latentGenotype_.at(locus);
    };

private:
    GenotypeMap<Data> observedGenotype_{};
    GenotypeMap<Parameter> latentGenotype_{};
};

#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
