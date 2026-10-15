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

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    class RandomAllelesBitSetSampler3 : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler3(std::shared_ptr<InfectionEventImpl> infection, std::shared_ptr<OrderingImpl> ordering, std::shared_ptr<InfectionEventImpl> latent_parent, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double fn_rate = 0.01, double fp_rate = 0.01) noexcept;

        void update() noexcept override;

        [[nodiscard]] Likelihood calculateSamplingProb() noexcept;

    private:
        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<OrderingImpl> ordering_;
        std::shared_ptr<InfectionEventImpl> latent_parent_;
        std::shared_ptr<T> target_;
        double fn_rate_;
        double fp_rate_;
        std::shared_ptr<Engine> rng_;

        ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize> proposalProbability_{fn_rate_, fp_rate_};
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::RandomAllelesBitSetSampler3(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<OrderingImpl> ordering,
            std::shared_ptr<InfectionEventImpl> latent_parent,
            std::shared_ptr<T> target,
            std::shared_ptr<Engine> rng,
            const double fn_rate,
            const double fp_rate) noexcept : infection_(std::move(infection)),
                                             ordering_(std::move(ordering)), latent_parent_(std::move(latent_parent)),
                                             target_(std::move(target)),
                                             fn_rate_(fn_rate),
                                             fp_rate_(fp_rate),
                                             rng_(rng) {}

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    void RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::update() noexcept {
        auto proposal = this->propose(*target_, *rng_);
        // before doing anything, calculate probability of sampling current state for adjustment
        const Likelihood current_state_prob = calculateSamplingProb();
        proposeGenotypesFromParentSet<MaxParentSetSize>(infection_, latent_parent_, ordering_->candidateParents(infection_), *rng_, fn_rate_, fp_rate_, [&](auto& genotype, const auto& value) { proposal.change(genotype, value); });
        const Likelihood proposed_state_prob = calculateSamplingProb();
        proposal.metropolis(current_state_prob - proposed_state_prob);
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    Likelihood RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, GeneticsImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::calculateSamplingProb() noexcept {
        // The probability of proposing the current genotypes (see ParentSetGenotypeProposal.h): a uniform mixture
        // over the parent-set options.
        return proposalProbability_.logProb(*infection_, *latent_parent_, ordering_->candidateParents(infection_));
    }

}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H
