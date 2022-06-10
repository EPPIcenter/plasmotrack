//
// Created by Maxwell Murphy on 1/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFECTION_H
#define TRANSMISSION_NETWORKS_APP_INFECTION_H


#include "core/abstract/observables/CheckpointablePassthrough.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"

#include "core/containers/Locus.h"

#include "core/datatypes/Data.h"

#include "core/parameters/Parameter.h"

#include <boost/container/flat_map.hpp>

#include <memory>
#include <utility>


namespace transmission_nets::core::containers {
    template<typename GeneticImpl, typename LocusImpl = Locus>
    class Infection : public abstract::Observable<Infection<GeneticImpl, LocusImpl>>,
                      public abstract::UncacheablePassthrough<Infection<GeneticImpl, LocusImpl>>,
                      public abstract::CheckpointablePassthrough<Infection<GeneticImpl, LocusImpl>> {


    public:
        template<typename Element>
        using GenotypeMap          = boost::container::flat_map<std::shared_ptr<LocusImpl>, Element>;// Maps genetic readout to loci
        using GenotypeDataMap      = GenotypeMap<std::shared_ptr<datatypes::Data<GeneticImpl>>>;
        using GenotypeParameterMap = GenotypeMap<std::shared_ptr<parameters::Parameter<GeneticImpl>>>;

        explicit Infection(std::string id, double observationTime);


        template<typename T>
        void addGenetics(std::shared_ptr<LocusImpl> locus, const T& obs, const T& latent);

        template<typename T>
        void addObservedGenetics(std::shared_ptr<LocusImpl> locus, const T& obs);

        template<typename T>
        void addLatentGenetics(std::shared_ptr<LocusImpl> locus, const T& latent);

        template<typename T>
        void addSourceGenetics(std::shared_ptr<LocusImpl> locus, const T& latent);

        GenotypeDataMap& observedGenotype() {
            return observedGenotype_;
        };

        GenotypeDataMap& observedGenotype() const {
            return observedGenotype_;
        };

        std::shared_ptr<datatypes::Data<GeneticImpl>> observedGenotype(std::shared_ptr<LocusImpl> locus) const {
            return observedGenotype_.at(locus);
        };

        std::shared_ptr<datatypes::Data<GeneticImpl>> observedGenotype(std::shared_ptr<LocusImpl> locus) {
            return observedGenotype_.at(locus);
        };

        GenotypeParameterMap& latentGenotype() {
            return latentGenotype_;
        };

        GenotypeParameterMap& latentGenotype() const {
            return latentGenotype_;
        };

        std::shared_ptr<parameters::Parameter<GeneticImpl>> latentGenotype(std::shared_ptr<LocusImpl> locus) {
            return latentGenotype_.at(locus);
        };

        std::shared_ptr<parameters::Parameter<GeneticImpl>> latentGenotype(std::shared_ptr<LocusImpl> locus) const {
            return latentGenotype_.at(locus);
        };

        GenotypeParameterMap& latentSource() {
            return latentSource_;
        };

        GenotypeParameterMap& latentSource() const {
            return latentSource_;
        };

        std::shared_ptr<parameters::Parameter<GeneticImpl>> latentSource(std::shared_ptr<LocusImpl> locus) {
            return latentSource_.at(locus);
        };

        std::shared_ptr<parameters::Parameter<GeneticImpl>> latentSource(std::shared_ptr<LocusImpl> locus) const {
            return latentSource_.at(locus);
        };

        const std::vector<std::shared_ptr<LocusImpl>>& loci() const {
            return loci_;
        }

        [[nodiscard]] std::shared_ptr<datatypes::Data<double>> observationTime() const {
            return observationTime_;
        }

        [[nodiscard]] std::shared_ptr<parameters::Parameter<double>> infectionDuration() const {
            return infectionDuration_;
        }

        double infectionTime() {
            double infectionTime = observationTime_->value() - infectionDuration_->value();
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
        GenotypeMap<std::shared_ptr<datatypes::Data<GeneticImpl>>> observedGenotype_{};
        GenotypeMap<std::shared_ptr<parameters::Parameter<GeneticImpl>>> latentGenotype_{};
        GenotypeMap<std::shared_ptr<parameters::Parameter<GeneticImpl>>> latentSource_{};
        std::vector<std::shared_ptr<LocusImpl>> loci_{};
        std::shared_ptr<datatypes::Data<double>> observationTime_;
        std::shared_ptr<parameters::Parameter<double>> infectionDuration_;// default to 100 days? or maybe something else -- look in constructor
    };

    template<typename GeneticImpl, typename LocusImpl>
    Infection<GeneticImpl, LocusImpl>::Infection(std::string id, const double observationTime) : id_(std::move(id)), observationTime_(std::make_shared<datatypes::Data<double>>(observationTime)) {
        infectionDuration_ = std::make_shared<parameters::Parameter<double>>(100);
        infectionDuration_->initializeValue(100.0);

        //        infectionDuration_->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        //        infectionDuration_->add_post_change_listener([=, this]() { this->notify_post_change(); });
        //        infectionDuration_->add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
        //        infectionDuration_->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        //        infectionDuration_->add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); })
    }

    //    template<typename GeneticImpl, typename LocusImpl>
    //    template<typename LocusDataIter>
    //    Infection<GeneticImpl, LocusImpl>::Infection(const std::string& id, const double observationTime, const LocusDataIter& obs, const LocusDataIter& latent) : Infection<GeneticImpl, LocusImpl>::Infection(id, observationTime) {
    //        for (const auto& [locus, genetics] : obs) {
    //            assert(locus->totalAlleles() == genetics.totalAlleles());
    //            observedGenotype_.emplace(locus, genetics);
    //        }
    //
    //        for (const auto& [locus, genetics] : latent) {
    //            loci_.push_back(locus);
    //            latentGenotype_.emplace(locus, genetics);
    //            // Creating pass through of notifications
    //            latentGenotype_.at(locus).add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
    //            latentGenotype_.at(locus).add_post_change_listener([=, this]() { this-> template notify_post_change(); });
    //            latentGenotype_.at(locus).add_save_state_listener([=, this](std::string savedStateId) { this-> template notify_save_state(savedStateId); });
    //            latentGenotype_.at(locus).add_accept_state_listener([=, this]() { this-> template notify_accept_state(); });
    //            latentGenotype_.at(locus).add_restore_state_listener([=, this](std::string savedStateId) { this-> template notify_restore_state(savedStateId); });
    //        }
    //    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(std::shared_ptr<LocusImpl> locus, const T& obs, const T& latent) {
        loci_.push_back(locus);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(latent);
        auto ob  = std::make_shared<datatypes::Data<GeneticImpl>>(obs);
        observedGenotype_.template insert_or_assign(locus, ob);
        latentGenotype_.template insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->template notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](std::string savedStateId) { this->template notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addObservedGenetics(std::shared_ptr<LocusImpl> locus, const T& obs) {
        /*
         * Copy the observed genetic state into the latent genetic state and latent source
         */
        loci_.push_back(locus);
        auto ob     = std::make_shared<datatypes::Data<GeneticImpl>>(obs);
        auto lat    = std::make_shared<parameters::Parameter<GeneticImpl>>(obs);
        auto latSrc = std::make_shared<parameters::Parameter<GeneticImpl>>(obs);
        observedGenotype_.template insert_or_assign(locus, ob);
        latentGenotype_.template insert_or_assign(locus, lat);
        latentSource_.template insert_or_assign(locus, latSrc);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->template notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](std::string savedStateId) { this->template notify_restore_state(savedStateId); });

        latentSource_.at(locus)->add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
        latentSource_.at(locus)->add_post_change_listener([=, this]() { this->template notify_post_change(); });
        latentSource_.at(locus)->add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
        latentSource_.at(locus)->add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
        latentSource_.at(locus)->add_restore_state_listener([=, this](std::string savedStateId) { this->template notify_restore_state(savedStateId); });
    }

    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addLatentGenetics(std::shared_ptr<LocusImpl> locus, const T& latent) {
        loci_.push_back(locus);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(latent);
        latentGenotype_.template insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->template notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](std::string savedStateId) { this->template notify_restore_state(savedStateId); });
    }


    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addSourceGenetics(std::shared_ptr<LocusImpl> locus, const T& latent) {
        loci_.push_back(locus);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(latent);
        latentSource_.template insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentSource_.at(locus)->add_pre_change_listener([=, this]() { this->template notify_pre_change(); });
        latentSource_.at(locus)->add_post_change_listener([=, this]() { this->template notify_post_change(); });
        latentSource_.at(locus)->add_save_state_listener([=, this](std::string savedStateId) { this->template notify_save_state(savedStateId); });
        latentSource_.at(locus)->add_accept_state_listener([=, this]() { this->template notify_accept_state(); });
        latentSource_.at(locus)->add_restore_state_listener([=, this](std::string savedStateId) { this->template notify_restore_state(savedStateId); });
    }

}// namespace transmission_nets::core::containers


#endif//TRANSMISSION_NETWORKS_APP_INFECTION_H
