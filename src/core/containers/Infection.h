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
    Infection(const std::string id, const LocusDataIter obs, const LocusDataIter latent);

    Infection(const std::string id);

    template<typename T>
    void addGenetics(LocusImpl *locus, const T &obs, const T &latent);

    template<typename T>
    void addGenetics(LocusImpl *locus, const T &obs);

    GenotypeMap<Data> &observedGenotype() {
        return observedGenotype_;
    };

    GenotypeMap<Data> &observedGenotype() const {
        return observedGenotype_;
    };


    const Data<GeneticImpl> &observedGenotype(LocusImpl *locus) const {
        return observedGenotype_.at(locus);
    };

    const GenotypeMap<Parameter> &latentGenotype() const {
        return latentGenotype_;
    };

    Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) {
        return latentGenotype_.at(locus);
    };

    const Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) const {
        return latentGenotype_.at(locus);
    };

    const std::vector<LocusImpl *> &loci() const {
        return loci_;
    }

    const std::string id() const {
        return id_;
    }

    const std::string serialize() const {
        return id_;
    }

private:
    std::string id_;
    GenotypeMap<Data> observedGenotype_{};
    GenotypeMap<Parameter> latentGenotype_{};
    std::vector<LocusImpl *> loci_{};
};

template<typename GeneticImpl, typename LocusImpl>
template<typename LocusDataIter>
Infection<GeneticImpl, LocusImpl>::Infection(const std::string id, const LocusDataIter obs, const LocusDataIter latent) : id_(id) {
    for (const auto& [locus, genetics] : obs) {
        assert(locus->totalAlleles() == genetics.totalAlleles());
        observedGenotype_.emplace(locus, genetics);
    }

    for (const auto& [locus, genetics] : latent) {
        loci_.push_back(locus);
        latentGenotype_.emplace(locus, genetics);
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=]() { this->notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=]() { this->notify_save_state(); });
        latentGenotype_.at(locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
    }
}

template<typename GeneticImpl, typename LocusImpl>
template<typename T>
void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs, const T &latent) {
    loci_.push_back(locus);
    latentGenotype_.emplace(locus, GeneticImpl(latent));
    observedGenotype_.emplace(locus, GeneticImpl(obs));
    // Creating pass through of notifications
    latentGenotype_.at(locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
    latentGenotype_.at(locus).add_post_change_listener([=]() { this->notify_post_change(); });
    latentGenotype_.at(locus).add_save_state_listener([=]() { this->notify_save_state(); });
    latentGenotype_.at(locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
    latentGenotype_.at(locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
}

template<typename GeneticImpl, typename LocusImpl>
Infection<GeneticImpl, LocusImpl>::Infection(const std::string id) : id_(id) {}

template<typename GeneticImpl, typename LocusImpl>
template<typename T>
void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs) {
    loci_.push_back(locus);
    latentGenotype_.emplace(locus, GeneticImpl(obs));
    observedGenotype_.emplace(locus, GeneticImpl(obs));
    latentGenotype_.at(locus).add_pre_change_listener([=]() { this->notify_pre_change(); });
    latentGenotype_.at(locus).add_post_change_listener([=]() { this->notify_post_change(); });
    latentGenotype_.at(locus).add_save_state_listener([=]() { this->notify_save_state(); });
    latentGenotype_.at(locus).add_accept_state_listener([=]() { this->notify_accept_state(); });
    latentGenotype_.at(locus).add_restore_state_listener([=]() { this->notify_restore_state(); });
}

#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
