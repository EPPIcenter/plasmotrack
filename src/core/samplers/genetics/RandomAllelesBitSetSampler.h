//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H

#include <boost/random.hpp>

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

        [[nodiscard]] float acceptanceRate() noexcept;

        [[nodiscard]] Likelihood logMetropolisHastingsAdjustment(const AllelesBitSetImpl& curr, const AllelesBitSetImpl& prop) noexcept;

    private:
        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int max_coi_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> allele_index_sampling_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
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
        SAMPLER_STATE_ID stateId = RandomAlleleBitSetID;
        Likelihood curLik = target_->value();
        parameter_->saveState(stateId);

        const auto& currentVal = parameter_->value();
        const auto& proposal = sampleProposal(currentVal);

        std::string currentValStr = currentVal.serialize();
        std::string proposalStr = proposal.serialize();

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const Likelihood adj = logMetropolisHastingsAdjustment(currentVal, proposal);
        const Likelihood acceptanceRatio = target_->value() - curLik + adj;
        const Likelihood logProbAccept = log(uniform_dist_(*rng_));
        const bool accept = (logProbAccept <= acceptanceRatio) and std::abs(acceptanceRatio) > 1e-10;


        if (debug_) {
            fmt::print("Acceptance ratio: {}\n", acceptanceRatio);
            fmt::print("Llik change: {}\n", std::abs(acceptanceRatio) > 1e-10);
            fmt::print("Log prob accept: {}\n", logProbAccept);
            fmt::print("Accept: {}\n", accept);
            fmt::print("Current likelihood: {}\n", curLik);
            fmt::print("Proposed likelihood: {}\n", target_->value());
            fmt::print("Adjustment: {}\n", adj);
            fmt::print("Current value: {}\n", currentValStr);
            fmt::print("Proposal: {}\n", proposalStr);
        }


        if (accept) {
            if (target_->value() < curLik and debug_) {
                fmt::print("Likelihood decreased!\n");
                fmt::print("Log prob accept: {}\n", logProbAccept);
                fmt::print("Current Likelihood: {}\n", curLik);
                fmt::print("New Llik: {} -- (AccR: {}, Adj: {})\n", target_->value(), acceptanceRatio, adj);
                fmt::print("Current value: {}\n", currentVal.serialize());
                fmt::print("Proposal: {}\n", proposal.serialize());

                parameter_->setValue(currentVal);
                fmt::print("Old Llik: {}\n", target_->value());
                parameter_->setValue(proposal);
                fmt::print("New Llik: {}\n", target_->value());
                fmt::print("------------------------\n");
            }
            acceptances_++;
            parameter_->acceptState();
        } else {
            rejections_++;
            parameter_->restoreState(stateId);
            assert(parameter_->value() == currentVal);
        }

        if (debug_) {
            fmt::print("Acceptances: {}\n", acceptances_);
            fmt::print("Rejections: {}\n", rejections_);
            fmt::print("------------------------\n");
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
    float RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptanceRate() noexcept {
        return float(acceptances_) / (acceptances_ + rejections_);
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
    Likelihood RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::logMetropolisHastingsAdjustment([[maybe_unused]] const AllelesBitSetImpl& curr, [[maybe_unused]] const AllelesBitSetImpl& prop) noexcept {
        // if there is only 1 allele, then there are total_alleles - 1 possible states to transition to.
        // If there are > 1 allele, then there are total_alleles possible states.
        // log(curr|proposed) - log(proposed|curr)
        unsigned int curr_total_alleles = curr.totalPositiveCount();
        unsigned int prop_total_alleles = prop.totalPositiveCount();

        Likelihood numerator = prop_total_alleles == 1 ? std::log(1) - std::log(curr.totalAlleles() - 1) : std::log(1) - std::log(curr.totalAlleles());
        Likelihood denominator = curr_total_alleles == 1 ? std::log(1) - std::log(prop.totalAlleles() - 1) : std::log(1) - std::log(prop.totalAlleles());

        return numerator - denominator;
        //        return 0;
    }
}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
