//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER2_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER2_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

#include <memory>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    class RandomAllelesBitSetSampler2 : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler2(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<LocusImpl> locus, std::shared_ptr<ParentSetImpl> ps, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int max_coi) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(AllelesBitSetImpl curr) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

        [[nodiscard]] Likelihood logMetropolisHastingsAdjustment(const AllelesBitSetImpl& curr, const AllelesBitSetImpl& prop) noexcept;

    private:
        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<LocusImpl> locus_;
        std::shared_ptr<ParentSetImpl> ps_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int max_coi_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> allele_index_sampling_dist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::RandomAllelesBitSetSampler2(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<LocusImpl> locus, std::shared_ptr<ParentSetImpl> ps, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int max_coi) noexcept : parameter_(std::move(parameter)), locus_(std::move(locus)), ps_(std::move(ps)), target_(std::move(target)), rng_(rng) , max_coi_(max_coi) {
        allele_index_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type(0, parameter_->value().totalAlleles() - 1));
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    void RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::update() noexcept {
        const std::string stateId = "State1";
        Likelihood curLik         = target_->value();
        parameter_->saveState(stateId);

        const auto currentVal = parameter_->value();
        const auto proposal = sampleProposal(parameter_->value());

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const auto adj = logMetropolisHastingsAdjustment(currentVal, proposal);
        const auto acceptanceRatio = target_->value() - curLik + adj;
        const auto logProbAccept   = log(uniform_dist_(*rng_));
        const bool accept          = logProbAccept <= acceptanceRatio;


//        fmt::print("Current Likelihood: {}\n", curLik);
//        fmt::print("New Llik: {} -- (AccR: {}, accepted: {})\n\n", target_->value(), acceptanceRatio, accept);

        if (accept) {
            acceptances_++;
            parameter_->acceptState();
        } else {
            rejections_++;
            parameter_->restoreState(stateId);
        }

        assert(!target_->isDirty());

        total_updates_++;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    unsigned int RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    unsigned int RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    double RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    AllelesBitSetImpl
    RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::sampleProposal(const AllelesBitSetImpl curr) noexcept {

        std::vector<unsigned int> valid_allele_indices{};
        auto tmp = curr;

        // Add alleles that can be flipped off
        for (size_t i = 0; i < curr.totalAlleles(); i++) {
            if (curr.allele(i)) {
                valid_allele_indices.push_back(i);
            }
        }


        // Add alleles that can be flipped on
        for (const auto& p : ps_->value()) {
            const auto latentGenotype = p->latentGenotype(locus_)->value();
            for (size_t i = 0; i < latentGenotype.totalAlleles(); i++) {
                if (latentGenotype.allele(i) && !curr.allele(i)) {
                    valid_allele_indices.push_back(i);
                }
            }
        }

        if (valid_allele_indices.size() <= 1) {
            return tmp;
        }

        bool flipped = false;
        while ((tmp.totalPositiveCount() == 0 or tmp.totalPositiveCount() > max_coi_) and !flipped) {
            flipped = false;
            int idx = allele_index_sampling_dist_(*rng_);
            if (std::find(valid_allele_indices.begin(), valid_allele_indices.end(), idx) != valid_allele_indices.end()) {
                tmp.flip(idx);
                flipped = true;
            }
        }

        return tmp;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, typename LocusImpl, typename ParentSetImpl>
    Likelihood RandomAllelesBitSetSampler2<T, Engine, AllelesBitSetImpl, LocusImpl, ParentSetImpl>::logMetropolisHastingsAdjustment(const AllelesBitSetImpl& curr, const AllelesBitSetImpl& prop) noexcept {
        // todo: fix the upper bound constraint also
        // if there is only 1 allele, then there are total_alleles - 1 possible states to transition to.
        // If there are > 1 allele, then there are total_alleles possible states.
        // log(curr|proposed) - log(proposed|curr)
        unsigned int curr_total_alleles = curr.totalPositiveCount();
        unsigned int prop_total_alleles = prop.totalPositiveCount();

        Likelihood numerator = prop_total_alleles == 1 ? std::log(1) - std::log(curr.totalAlleles() - 1) : std::log(1) - std::log(curr.totalAlleles());
        Likelihood denominator = curr_total_alleles == 1 ? std::log(1) - std::log(prop.totalAlleles() - 1) : std::log(1) - std::log(prop.totalAlleles());

        return numerator - denominator;
    }
}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER2_H
