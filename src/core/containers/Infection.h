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
        Infection(const std::string& id, double observationTime, LocusDataIter obs, LocusDataIter latent);

        explicit Infection(const std::string& id, double observationTime);

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

       datatypes::Data<double> &observationTime() {
           return observationTime_;
       }

       parameters::Parameter<double> &infectionDuration() {
           return infectionDuration_;
       }

       double infectionTime() {
           double infectionTime = observationTime_.value() - infectionDuration_.value();
           return infectionTime;
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
        datatypes::Data<double> observationTime_;
        parameters::Parameter<double> infectionDuration_; // default to 10 days? or maybe something else -- look in constructor

    };

    template<typename GeneticImpl, typename LocusImpl>
    Infection<GeneticImpl, LocusImpl>::Infection(const std::string& id, const double observationTime) : id_(id), observationTime_(observationTime) {
        // pass through of notifications to listeners of "infection"
        infectionDuration_.initializeValue(100.0);
//        infectionDuration_.add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
//        infectionDuration_.add_post_change_listener([=, this]() { this->template notify_post_change(); });
//        infectionDuration_.add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
//        infectionDuration_.add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
//        infectionDuration_.add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename LocusDataIter>
    Infection<GeneticImpl, LocusImpl>::Infection(const std::string& id, const double observationTime, const LocusDataIter obs, const LocusDataIter latent) : Infection<GeneticImpl, LocusImpl>::Infection(id, observationTime) {
        for (const auto& [locus, genetics] : obs) {
            assert(locus->totalAlleles() == genetics.totalAlleles());
            observedGenotype_.emplace(locus, genetics);
        }

        for (const auto& [locus, genetics] : latent) {
            loci_.push_back(locus);
            latentGenotype_.emplace(locus, genetics);
            // Creating pass through of notifications
            latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
            latentGenotype_.at(locus).add_post_change_listener([=, this]() { this-> template notify_post_change(); });
            latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this-> template notify_save_state(savedStateId); });
            latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this-> template notify_accept_state(); });
            latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
        }
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs, const T &latent) {
        loci_.push_back(locus);
        observedGenotype_.emplace(locus, GeneticImpl(obs));
        latentGenotype_.emplace(locus, GeneticImpl(latent));
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this-> template notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this-> template notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this-> template notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this-> template notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(LocusImpl *locus, const T &obs) {
        loci_.push_back(locus);
        observedGenotype_.emplace(locus, GeneticImpl(obs));
        latentGenotype_.emplace(locus, GeneticImpl(obs));
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this-> template notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this-> template notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this-> template notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this-> template notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addLatentGenetics(LocusImpl *locus, const T &latent) {
        loci_.push_back(locus);
        latentGenotype_.emplace(locus, GeneticImpl(latent));
        // Creating pass through of notifications
        latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this-> template notify_pre_change(); });
        latentGenotype_.at(locus).add_post_change_listener([=, this]() { this-> template notify_post_change(); });
        latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this-> template notify_save_state(savedStateId); });
        latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this-> template notify_accept_state(); });
        latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
    }
}



#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
