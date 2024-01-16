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

        explicit Infection(std::string id, float observationTime, bool symptomatic = true);

        Infection(const Infection& other, const std::string& id = "", bool retain_alleles = true) {
            static int uid = 0;
            uid_ = uid++;
            // Copy constructor -- create a new infection from an existing one using the latent genetics.
            if (id.empty()) {
                id_ = other.id_ + "_copy";
            } else {
                id_ = id;
            }

            observationTime_   = other.observationTime_;
            infectionDuration_ = other.infectionDuration_;
            symptomatic_       = other.symptomatic_;

            for (const auto& [locus, data] : other.latentGenotype_) {
                if (retain_alleles) {
                    addGenetics(locus, data->value(), data->value());
                } else {
                    const auto obsGenetics = data->value().allelesStr();
                    std::string bitstr;
                    int alleleFound = 0;
                    for (const auto c : obsGenetics) {
                        if (c == '1' and alleleFound < 2) {
                            bitstr += "1";
                            alleleFound += true;
                        } else {
                            bitstr += "0";
                        }
                    }
                    if (!alleleFound) {
                        bitstr[bitstr.size() - 1] = '1';
                    }
                    addGenetics(locus, GeneticImpl(bitstr), GeneticImpl(bitstr));
                }
            }
        }


        template<typename T>
        void addGenetics(std::shared_ptr<LocusImpl> locus, const T& obs, const T& latent);

        template<typename T>
        void addObservedGenetics(std::shared_ptr<LocusImpl> locus, const T& obs);

        template<typename T>
        void addLatentGenetics(std::shared_ptr<LocusImpl> locus, const T& latent);

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

        /**
         * @brief Returns a vector of all the loci in the infection.
         * @return A vector of all the loci in the infection.
         */
        const std::vector<std::shared_ptr<LocusImpl>>& loci() const {
            return loci_;
        }

        /**
         * @brief Returns the observation time of the infection.
         * @return The observation time of the infection.
         */
        [[nodiscard]] std::shared_ptr<datatypes::Data<float>> observationTime() const {
            return observationTime_;
        }

        /**
         * @brief Returns the duration of the infection.
         * @return The duration of the infection.
         */
        [[nodiscard]] std::shared_ptr<parameters::Parameter<float>> infectionDuration() const {
            return infectionDuration_;
        }

        /**
         * @brief Returns the time of the infection.
         * @return The time of the infection.
         */
        float infectionTime() {
            float infectionTime = observationTime_->value() - infectionDuration_->value();
            return infectionTime;
        }

        /**
         * @brief Sets the id of the infection.
         * @param id The id of the infection.
         */
        void setId(const std::string& id) {
            id_ = id;
        }

        /**
         * @brief Returns the id of the infection.
         * @return The id of the infection.
         */
        [[nodiscard]] std::string id() const {
            return id_;
        }

        [[nodiscard]] int uid() const {
            return uid_;
        }

        /**
         * @brief Returns the string representation of the infection using the id.
         * @return The string representation of the infection.
         */
        [[nodiscard]] std::string serialize() const {
            return id_;
        }


        // Infections are comparable by their id.
        bool operator<(const Infection& rhs) const {
            if (static_cast<const transmission_nets::core::abstract::Observable<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this) < static_cast<const transmission_nets::core::abstract::Observable<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs))
                return true;
            if (static_cast<const transmission_nets::core::abstract::Observable<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs) < static_cast<const transmission_nets::core::abstract::Observable<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this))
                return false;
            if (static_cast<const transmission_nets::core::abstract::UncacheablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this) < static_cast<const transmission_nets::core::abstract::UncacheablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs))
                return true;
            if (static_cast<const transmission_nets::core::abstract::UncacheablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs) < static_cast<const transmission_nets::core::abstract::UncacheablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this))
                return false;
            if (static_cast<const transmission_nets::core::abstract::CheckpointablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this) < static_cast<const transmission_nets::core::abstract::CheckpointablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs))
                return true;
            if (static_cast<const transmission_nets::core::abstract::CheckpointablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(rhs) < static_cast<const transmission_nets::core::abstract::CheckpointablePassthrough<transmission_nets::core::containers::Infection<GeneticImpl, LocusImpl>>&>(*this))
                return false;
            return id_ < rhs.id_;
        }

        bool operator>(const Infection& rhs) const {
            return rhs < *this;
        }

        bool operator<=(const Infection& rhs) const {
            return !(rhs < *this);
        }

        bool operator>=(const Infection& rhs) const {
            return !(*this < rhs);
        }

        bool isSymptomatic() const {
            return symptomatic_->value();
        }

    private:
        std::string id_;
        int uid_;
        GenotypeMap<std::shared_ptr<datatypes::Data<GeneticImpl>>> observedGenotype_{};
        GenotypeMap<std::shared_ptr<parameters::Parameter<GeneticImpl>>> latentGenotype_{};
        std::vector<std::shared_ptr<LocusImpl>> loci_{};
        std::shared_ptr<datatypes::Data<float>> observationTime_;
        std::shared_ptr<parameters::Parameter<float>> infectionDuration_;// default to 100 days? or maybe something else -- look in constructor
        std::shared_ptr<datatypes::Data<bool>> symptomatic_{};
    };

    template<typename GeneticImpl, typename LocusImpl>
    Infection<GeneticImpl, LocusImpl>::Infection(std::string id, const float observationTime, const bool symptomatic) : id_(std::move(id)), observationTime_(std::make_shared<datatypes::Data<float>>(observationTime)), symptomatic_(std::make_shared<datatypes::Data<bool>>(symptomatic)) {
        static int uid = 0;
        uid_ = uid++;
        infectionDuration_ = std::make_shared<parameters::Parameter<float>>(10.0);
        infectionDuration_->initializeValue(10.0);

        //        infectionDuration_->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        //        infectionDuration_->add_post_change_listener([=, this]() { this->notify_post_change(); });
        //        infectionDuration_->add_save_state_listener([=, this](std::string savedStateId) { this->notify_save_state(savedStateId); });
        //        infectionDuration_->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        //        infectionDuration_->add_restore_state_listener([=, this](std::string savedStateId) { this->notify_restore_state(savedStateId); })
    }



    /**
     * @brief Add observed and latent genotypes to the infection.
     * @tparam GeneticImpl Class implementing the Genetic interface.
     * @tparam LocusImpl Class implementing the Locus interface.
     * @tparam T Class that may be used to initialize the genotypes.
     * @param locus Pointer to the locus.
     * @param obs The observed genotype.
     * @param latent The latent genotype.
     */
    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addGenetics(std::shared_ptr<LocusImpl> locus, const T& obs, const T& latent) {
        loci_.push_back(locus);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(latent);
        auto ob  = std::make_shared<datatypes::Data<GeneticImpl>>(obs);
        observedGenotype_.insert_or_assign(locus, ob);
        latentGenotype_.insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](int savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](int savedStateId) { this->notify_restore_state(savedStateId); });
    }

    /**
     * @brief Add observed and latent genotypes to the infection. Latent genotypes are initialized to the observed genotypes.
     * @tparam GeneticImpl Class implementing the Genetic interface.
     * @tparam LocusImpl Class implementing the Locus interface.
     * @tparam T Class that may be used to initialize the genotypes.
     * @param locus Pointer to the locus.
     * @param obs The observed genotype.
     */
    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addObservedGenetics(std::shared_ptr<LocusImpl> locus, const T& obs) {
        loci_.push_back(locus);
        auto ob  = std::make_shared<datatypes::Data<GeneticImpl>>(obs);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(obs);
        observedGenotype_.insert_or_assign(locus, ob);
        latentGenotype_.insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](int savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](int savedStateId) { this->notify_restore_state(savedStateId); });
    }


    /**
     * @brief Add latent genotypes to the infection.
     * @tparam GeneticImpl Class implementing the Genetic interface.
     * @tparam LocusImpl Class implementing the Locus interface.
     * @tparam T Class that may be used to initialize the genotypes.
     * @param locus Pointer to the locus.
     * @param latent The latent genotype.
     */
    template<typename GeneticImpl, typename LocusImpl>
    template<typename T>
    void Infection<GeneticImpl, LocusImpl>::addLatentGenetics(std::shared_ptr<LocusImpl> locus, const T& latent) {
        loci_.push_back(locus);
        auto lat = std::make_shared<parameters::Parameter<GeneticImpl>>(latent);
        latentGenotype_.insert_or_assign(locus, lat);
        // Creating pass through of notifications
        latentGenotype_.at(locus)->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        latentGenotype_.at(locus)->add_post_change_listener([=, this]() { this->notify_post_change(); });
        latentGenotype_.at(locus)->add_save_state_listener([=, this](int savedStateId) { this->notify_save_state(savedStateId); });
        latentGenotype_.at(locus)->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        latentGenotype_.at(locus)->add_restore_state_listener([=, this](int savedStateId) { this->notify_restore_state(savedStateId); });
    }

}// namespace transmission_nets::core::containers


#endif//TRANSMISSION_NETWORKS_APP_INFECTION_H
