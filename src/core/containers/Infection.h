//
// Created by Maxwell Murphy on 1/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFECTION_H
#define TRANSMISSION_NETWORKS_APP_INFECTION_H

#include <boost/container/flat_map.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"

#include "core/containers/Locus.h"

#include "core/datatypes/Data.h"

#include "core/parameters/Parameter.h"


template<typename GeneticImpl, typename LocusImpl = Locus>
class Infection : public Observable<Infection<GeneticImpl, LocusImpl>>,
                  public UncacheablePassthrough<Infection<GeneticImpl, LocusImpl>>,
                  public CheckpointablePassthrough<Infection<GeneticImpl, LocusImpl>> {

    template<template<typename> typename Wrapper>
    using GenotypeMap = boost::container::flat_map<LocusImpl *, Wrapper<GeneticImpl>>; // Maps genetic readout to loci

public:
    using LocusGeneticsAssignment = std::pair<LocusImpl *, GeneticImpl>;

    template<typename LocusDataIter>
    Infection(LocusDataIter obs, LocusDataIter latent);

    Infection();

    template<typename T>
    void addGenetics(LocusImpl *locus, const T &obs, const T &latent);

    GenotypeMap<Data> &observedGenotype() {
        return observedGenotype_;
    };

    GenotypeMap<Data> &observedGenotype() const {
        return observedGenotype_;
    };

    Data<GeneticImpl> &observedGenotype(LocusImpl *locus) {
        return observedGenotype_.at(locus);
    };

    Data<GeneticImpl> &observedGenotype(LocusImpl *locus) const {
        return observedGenotype_.at(locus);
    };

    GenotypeMap<Parameter> &latentGenotype() {
        return latentGenotype_;
    };

    GenotypeMap<Parameter> &latentGenotype() const {
        return latentGenotype_;
    };

    Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) {
        return latentGenotype_.at(locus);
    };

    Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) const {
        return latentGenotype_.at(locus);
    };

private:
    GenotypeMap<Data> observedGenotype_{};
    GenotypeMap<Parameter> latentGenotype_{};
};

template<typename GeneticImpl, typename LocusImpl>
template<typename LocusDataIter>
Infection<GeneticImpl, LocusImpl>::Infection(LocusDataIter obs, LocusDataIter latent) {
    for (auto const&[locus, genetics] : obs) {
        assert(locus->totalAlleles() == genetics.totalAlleles());
        observedGenotype_.emplace(locus, genetics);
    }

    // TODO: Generalize this notification forwarding in collections of parameters
    for (auto const&[locus, genetics] : latent) {
        latentGenotype_.emplace(locus, genetics);
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([&]() { this->notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([&]() { this->notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([&]() { this->notify_save_state(); });
        latentGenotype_.at(locus).add_accept_state_listener([&]() { this->notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([&]() { this->notify_restore_state(); });
    }
}

template<typename GeneticImpl, typename LocusImpl>
template<typename T>
void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs, const T &latent) {
    latentGenotype_.emplace(locus, GeneticImpl(latent));
    observedGenotype_.emplace(locus, GeneticImpl(obs));
    latentGenotype_.at(locus).add_pre_change_listener([&]() { this->notify_pre_change(); });
    latentGenotype_.at(locus).add_post_change_listener([&]() { this->notify_post_change(); });
    latentGenotype_.at(locus).add_save_state_listener([&]() { this->notify_save_state(); });
    latentGenotype_.at(locus).add_accept_state_listener([&]() { this->notify_accept_state(); });
    latentGenotype_.at(locus).add_restore_state_listener([&]() { this->notify_restore_state(); });
}

template<typename GeneticImpl, typename LocusImpl>
Infection<GeneticImpl, LocusImpl>::Infection() {}

#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
