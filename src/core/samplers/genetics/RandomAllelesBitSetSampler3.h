//
// Created by Maxwell Murphy on 10/18/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H

#include <boost/random.hpp>

#include "core/datatypes/Genotype.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/genetics/ObservationRate.h"
#include "core/samplers/genetics/ParentSetGenotypeProposal.h"
#include "core/samplers/genetics/ParentSetProposalProbability.h"
#include "core/utils/numerics.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    // An infection's parent-set log-likelihoods in parent-set-distribution order (see parentSetIndexOfOption()) at
    // the current state, written to the vector; returns their log-sum.
    using ParentSetLliks = std::function<double(std::vector<double>&)>;

    template<typename T, typename Engine, typename InfectionEventImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    class RandomAllelesBitSetSampler3 : public AbstractSampler {
    public:
        // `expectedFalsePositives` and `expectedFalseNegatives`, when set, replace fn_rate and fp_rate at each update
        // with those counts divided by the infection's mean alleles per locus.
        //
        // With `parentSetLliks`, the parent-set option is drawn from a mixture: with probability `uniformShare`
        // uniformly, otherwise from the infection's parent-set posterior at the current state. The posterior changes
        // with the child's genotypes, so the reverse move's weights are read at the proposed state. With
        // `lociPerUpdate` above 0, each update redraws only that many loci, chosen uniformly; the rest are kept.
        RandomAllelesBitSetSampler3(std::shared_ptr<InfectionEventImpl> infection, std::shared_ptr<OrderingImpl> ordering, std::shared_ptr<InfectionEventImpl> latent_parent, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double fn_rate = 0.01, double fp_rate = 0.01, std::shared_ptr<core::parameters::Parameter<double>> expectedFalsePositives = nullptr, std::shared_ptr<core::parameters::Parameter<double>> expectedFalseNegatives = nullptr, ParentSetLliks parentSetLliks = nullptr, double uniformShare = 1.0, size_t lociPerUpdate = 0) noexcept;

        void update() noexcept override;

        [[nodiscard]] Likelihood calculateSamplingProb() noexcept;

    private:
        // The option probabilities at the current state, in parent-set-distribution order.
        void optionProbs(std::vector<double>& out);

        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<OrderingImpl> ordering_;
        std::shared_ptr<InfectionEventImpl> latent_parent_;
        std::shared_ptr<T> target_;
        double fn_rate_;
        double fp_rate_;
        std::shared_ptr<Engine> rng_;
        std::shared_ptr<core::parameters::Parameter<double>> expectedFp_;
        std::shared_ptr<core::parameters::Parameter<double>> expectedFn_;
        ParentSetLliks parentSetLliks_;
        double uniformShare_;
        size_t lociPerUpdate_;
        std::vector<const void*> chosenLoci_{};
        std::vector<size_t> lociOrder_{};

        ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize> proposalProbability_{fn_rate_, fp_rate_};
        typename ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::Evaluation currentEvaluation_{};
        typename ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::Evaluation proposedEvaluation_{};
        std::vector<double> currentProbs_{};
        std::vector<double> proposedProbs_{};
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::RandomAllelesBitSetSampler3(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<OrderingImpl> ordering,
            std::shared_ptr<InfectionEventImpl> latent_parent,
            std::shared_ptr<T> target,
            std::shared_ptr<Engine> rng,
            const double fn_rate,
            const double fp_rate,
            std::shared_ptr<core::parameters::Parameter<double>> expectedFalsePositives,
            std::shared_ptr<core::parameters::Parameter<double>> expectedFalseNegatives,
            ParentSetLliks parentSetLliks,
            const double uniformShare,
            const size_t lociPerUpdate) noexcept : infection_(std::move(infection)),
                                             ordering_(std::move(ordering)), latent_parent_(std::move(latent_parent)),
                                             target_(std::move(target)),
                                             fn_rate_(fn_rate),
                                             fp_rate_(fp_rate),
                                             rng_(rng),
                                             expectedFp_(std::move(expectedFalsePositives)),
                                             expectedFn_(std::move(expectedFalseNegatives)),
                                             parentSetLliks_(std::move(parentSetLliks)),
                                             uniformShare_(uniformShare),
                                             lociPerUpdate_(lociPerUpdate) {}

    template<typename T, typename Engine, typename InfectionEventImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    void RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::update() noexcept {
        if (expectedFp_ and expectedFn_) {
            double alleles = 0;
            int loci       = 0;
            for (const auto& locus : infection_->loci()) {
                alleles += locus->totalAlleles();
                ++loci;
            }
            const double mean = loci == 0 ? 1.0 : alleles / static_cast<double>(loci);
            fn_rate_ = perAlleleRate(expectedFn_->value(), mean);
            fp_rate_ = perAlleleRate(expectedFp_->value(), mean);
            proposalProbability_.setRates(fn_rate_, fp_rate_);
        }
        const auto& parents = ordering_->candidateParents(infection_);
        auto change         = [&](auto& proposal) { return [&](auto& genotype, const auto& value) { proposal.change(genotype, value); }; };

        const bool weighted = parentSetLliks_ and uniformShare_ < 1.0;
        const auto& loci    = infection_->loci();
        const bool subset   = lociPerUpdate_ > 0 and lociPerUpdate_ < loci.size();
        if (not weighted and not subset) {
            auto proposal = this->propose(*target_, *rng_);
            // before doing anything, calculate probability of sampling current state for adjustment
            const Likelihood current_state_prob = calculateSamplingProb();
            proposeGenotypesFromParentSet<MaxParentSetSize>(infection_, latent_parent_, parents, *rng_, fn_rate_, fp_rate_, change(proposal));
            proposalProbability_.evaluateChildOnly(*infection_, proposedEvaluation_);
            const Likelihood proposed_state_prob = proposalProbability_.combine(proposedEvaluation_, {});
            proposal.metropolis(current_state_prob - proposed_state_prob);
            return;
        }

        chosenLoci_.clear();
        if (subset) {
            lociOrder_.resize(loci.size());
            std::iota(lociOrder_.begin(), lociOrder_.end(), 0);
            for (size_t k = 0; k < lociPerUpdate_; ++k) {
                boost::random::uniform_int_distribution<size_t> pick(k, loci.size() - 1);
                std::swap(lociOrder_[k], lociOrder_[pick(*rng_)]);
                chosenLoci_.push_back(loci[lociOrder_[k]].get());
            }
        }
        auto include = [&](const auto& locus) { return not subset or std::ranges::find(chosenLoci_, static_cast<const void*>(locus.get())) != chosenLoci_.end(); };

        auto proposal = this->propose(*target_, *rng_);
        proposalProbability_.evaluate(*infection_, *latent_parent_, parents, currentEvaluation_, include);
        if (weighted) {
            optionProbs(currentProbs_);
        } else {
            currentProbs_.clear();
        }

        const size_t numCombinations = numParentCombinations<MaxParentSetSize>(parents.size());
        size_t option                = 0;
        if (weighted) {
            boost::random::uniform_01<> unif;
            double u            = unif(*rng_);
            size_t parentSetIdx = currentProbs_.size() - 1;
            for (size_t i = 0; i < currentProbs_.size(); ++i) {
                if (u < currentProbs_[i]) {
                    parentSetIdx = i;
                    break;
                }
                u -= currentProbs_[i];
            }
            option = optionOfParentSetIndex(parentSetIdx, numCombinations);
        } else {
            option = boost::random::uniform_int_distribution<size_t>(0, 2 * numCombinations)(*rng_);
        }
        proposeGenotypesFromOption<MaxParentSetSize>(infection_, latent_parent_, parents, option, *rng_, fn_rate_, fp_rate_, change(proposal), include);

        proposalProbability_.evaluateChildOnly(*infection_, proposedEvaluation_);
        if (weighted) {
            optionProbs(proposedProbs_);
        } else {
            proposedProbs_.clear();
        }
        const double reverse = proposalProbability_.combine(currentEvaluation_, proposedProbs_);
        const double forward = proposalProbability_.combine(proposedEvaluation_, currentProbs_);
        proposal.metropolis(reverse - forward);
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    void RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::optionProbs(std::vector<double>& out) {
        const double total = parentSetLliks_(out);
        const double share = uniformShare_ / static_cast<double>(out.size());
        // A state the posterior rules out is rejected whatever the weights, so uniform ones stand in.
        const bool uniform = not std::isfinite(total);
        for (auto& p : out) {
            p = uniform ? 1.0 / static_cast<double>(out.size()) : share + (1.0 - uniformShare_) * std::exp(p - total);
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename OrderingImpl, int MaxParentSetSize, int MaxCOI>
    Likelihood RandomAllelesBitSetSampler3<T, Engine, InfectionEventImpl, OrderingImpl, MaxParentSetSize, MaxCOI>::calculateSamplingProb() noexcept {
        // The probability of proposing the current genotypes (see ParentSetGenotypeProposal.h): a uniform mixture
        // over the parent-set options.
        return proposalProbability_.logProb(*infection_, *latent_parent_, ordering_->candidateParents(infection_));
    }

}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ALLELESAMPLER3_H
