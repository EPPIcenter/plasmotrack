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

namespace transmission_nets::core::containers {
    template<typename GeneticImpl, typename LocusImpl = Locus>
    class Infection : public abstract::Observable<Infection<GeneticImpl, LocusImpl>>,
                      public abstract::UncacheablePassthrough<Infection<GeneticImpl, LocusImpl>>,
                      public abstract::CheckpointablePassthrough<Infection<GeneticImpl, LocusImpl>> {

        template<template<typename> typename Wrapper>
        using GenotypeMap = boost::container::flat_map<LocusImpl *, Wrapper<GeneticImpl>>; // Maps genetic readout to loci

    public:
        using LocusGeneticsAssignment = std::pair<LocusImpl *, GeneticImpl>;

        template<typename LocusDataIter>
        Infection(std::string id, LocusDataIter obs, LocusDataIter latent);

        explicit Infection(std::string id);

        template<typename T>
        void addGenetics(LocusImpl *locus, const T &obs, const T &latent);

        template<typename T>
        void addGenetics(LocusImpl *locus, const T &obs);

        template<typename T>
        void addLatentGenetics(LocusImpl *locus, const T &latent);

        GenotypeMap<datatypes::Data> &observedGenotype() {
            return observedGenotype_;
        };

        GenotypeMap<datatypes::Data> &observedGenotype() const {
            return observedGenotype_;
        };


        const datatypes::Data<GeneticImpl> &observedGenotype(LocusImpl *locus) const {
            return observedGenotype_.at(locus);
        };

        const GenotypeMap<parameters::Parameter> &latentGenotype() const {
            return latentGenotype_;
        };

        parameters::Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) {
            return latentGenotype_.at(locus);
        };

        const parameters::Parameter<GeneticImpl> &latentGenotype(LocusImpl *locus) const {
            return latentGenotype_.at(locus);
        };

        const std::vector<LocusImpl *> &loci() const {
            return loci_;
        }

        [[nodiscard]] std::string id() const {
            return id_;
        }

        [[nodiscard]] std::string serialize() const {
            return id_;
        }

    private:
        std::string id_;
        GenotypeMap<datatypes::Data> observedGenotype_{};
        GenotypeMap<parameters::Parameter> latentGenotype_{};
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
            latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->notify_pre_change(); });
            latentGenotype_.at(locus).add_post_change_listener([=, this]() { this->notify_post_change(); });
            latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
            latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this->notify_accept_state(); });
            latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); });
        }
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs, const T &latent) {
        loci_.push_back(locus);
        observedGenotype_.emplace(locus, GeneticImpl(obs));
        latentGenotype_.emplace(locus, GeneticImpl(latent));
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    Infection<GeneticImpl, LocusImpl>::Infection(const std::string id) : id_(id) {}

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs) {
        loci_.push_back(locus);
        observedGenotype_.emplace(locus, GeneticImpl(obs));
        latentGenotype_.emplace(locus, GeneticImpl(obs));
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addLatentGenetics(LocusImpl *locus, const T &latent) {
        loci_.push_back(locus);
        latentGenotype_.emplace(locus, GeneticImpl(latent));
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); });
    }
}



#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
