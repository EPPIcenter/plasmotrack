//
// Created by Maxwell Murphy on 10/18/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H

#include <boost/random.hpp>

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/genetics/ParentSetGenotypeProposal.h"
#include "core/samplers/genetics/ParentSetProposalProbability.h"
#include "core/utils/numerics.h"

#include <memory>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    class RandomAllelesBitSetSampler3 : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler3(std::shared_ptr<InfectionEventImpl> infection, std::shared_ptr<ParentSetImpl> ps, std::shared_ptr<InfectionEventImpl> latent_parent, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double fn_rate = 0.01, double fp_rate = 0.01) noexcept;

        void update() noexcept override;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

        [[nodiscard]] Likelihood calculateSamplingProb() noexcept;

    private:
        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<ParentSetImpl> ps_;
        std::shared_ptr<InfectionEventImpl> latent_parent_;
        std::shared_ptr<T> target_;
        double fn_rate_;
        double fp_rate_;
        std::shared_ptr<Engine> rng_;

        boost::random::uniform_01<> uniform_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

        ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize> proposalProbability_{fn_rate_, fp_rate_};
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::RandomAllelesBitSetSampler3(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<ParentSetImpl> ps,
            std::shared_ptr<InfectionEventImpl> latent_parent,
            std::shared_ptr<T> target,
            std::shared_ptr<Engine> rng,
            const double fn_rate,
            const double fp_rate) noexcept : infection_(std::move(infection)),
                                             ps_(std::move(ps)), latent_parent_(std::move(latent_parent)),
                                             target_(std::move(target)),
                                             fn_rate_(fn_rate),
                                             fp_rate_(fp_rate),
                                             rng_(rng) {}

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    void RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::update() noexcept {
        SAMPLER_STATE_ID stateId = RandomAlleleBitSet3ID;

        Likelihood cur_lik = target_->value();
        // before doing anything, calculate probability of sampling current state for adjustment
        Likelihood current_state_prob = calculateSamplingProb();

        proposeGenotypesFromParentSet<MaxParentSetSize>(infection_, latent_parent_, ps_->value(), *rng_, fn_rate_, fp_rate_, stateId);

        Likelihood proposed_state_prob = calculateSamplingProb();

        const auto acceptanceRatio = target_->value() - cur_lik + current_state_prob - proposed_state_prob;
        const auto logProbAccept = log(uniform_dist_(*rng_));
        const bool accept = logProbAccept <= acceptanceRatio;

        if (accept) {
            acceptances_++;
            countProposal(true);
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->acceptState();
            }
        } else {
            rejections_++;
            countProposal(false);
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->restoreState(stateId);
            }
        }

        assert(!target_->isDirty());

        total_updates_++;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    Likelihood RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::calculateSamplingProb() noexcept {
        // The probability of proposing the current genotypes (see ParentSetGenotypeProposal.h): a uniform mixture
        // over the parent-set options.
        return proposalProbability_.logProb(*infection_, *latent_parent_, ps_->value());
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    unsigned int RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    unsigned int RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    double RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::acceptanceRate() noexcept {
        return static_cast<double>(acceptances_) / (acceptances_ + rejections_);
    }

}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H
