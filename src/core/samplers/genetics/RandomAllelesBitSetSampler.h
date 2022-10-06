//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

#include <memory>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    class RandomAllelesBitSetSampler : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int max_coi) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(AllelesBitSetImpl curr) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

        [[nodiscard]] Likelihood logMetropolisHastingsAdjustment(const AllelesBitSetImpl& curr, const AllelesBitSetImpl& prop) noexcept;

    private:
        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int max_coi_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> allele_index_sampling_dist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::RandomAllelesBitSetSampler(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int max_coi) noexcept : parameter_(parameter), target_(target), rng_(rng), max_coi_(max_coi) {
        allele_index_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type(0, parameter_->value().totalAlleles() - 1));
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    void RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::update() noexcept {
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

//        assert(!target_->isDirty());

        total_updates_++;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    double RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    AllelesBitSetImpl
    RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::sampleProposal(const AllelesBitSetImpl curr) noexcept {
        auto tmp = curr;
        tmp.flip(allele_index_sampling_dist_(*rng_));

        while (tmp.totalPositiveCount() == 0 or tmp.totalPositiveCount() > max_coi_) {
            tmp.flip(allele_index_sampling_dist_(*rng_));
        }

        return tmp;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    Likelihood RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::logMetropolisHastingsAdjustment(const AllelesBitSetImpl& curr, const AllelesBitSetImpl& prop) noexcept {
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


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
