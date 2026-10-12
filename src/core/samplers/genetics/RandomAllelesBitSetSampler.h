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

        AllelesBitSetImpl sampleProposal(const AllelesBitSetImpl& curr) noexcept;

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

        const auto currentVal = parameter_->value();
        const auto& proposal = sampleProposal(currentVal);

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const Likelihood adj = logMetropolisHastingsAdjustment(currentVal, proposal);
        const Likelihood acceptanceRatio = target_->value() - curLik + adj;
        const Likelihood logProbAccept = log(uniform_dist_(*rng_));
        const bool accept = (logProbAccept <= acceptanceRatio) and std::abs(acceptanceRatio) > 1e-10;

        if (accept) {
            acceptances_++;
            countProposal(true);
            parameter_->acceptState();
        } else {
            rejections_++;
            countProposal(false);
            parameter_->restoreState(stateId);
            assert(parameter_->value() == currentVal);
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
    RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::sampleProposal(const AllelesBitSetImpl& curr) noexcept {
        auto tmp = curr;
        tmp.flip(allele_index_sampling_dist_(*rng_));

        while (tmp.totalPositiveCount() == 0 or tmp.totalPositiveCount() > max_coi_) {
            tmp.flip(allele_index_sampling_dist_(*rng_));
        }

        return tmp;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    Likelihood RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::logMetropolisHastingsAdjustment([[maybe_unused]] const AllelesBitSetImpl& curr, [[maybe_unused]] const AllelesBitSetImpl& prop) noexcept {
        // The proposal flips uniformly chosen alleles until the genotype is valid again. Each flip is as likely as its
        // reverse, so every path is as likely as the same path reversed and the proposal is symmetric.
        return 0;
    }
}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
