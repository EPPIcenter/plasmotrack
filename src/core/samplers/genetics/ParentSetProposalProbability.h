//
// The proposal probability of RandomAllelesBitSetSampler3.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETPROPOSALPROBABILITY_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETPROPOSALPROBABILITY_H

#include "core/samplers/genetics/GenomePartition.h"
#include "core/samplers/genetics/ParentSetGenotypeProposal.h"
#include "core/utils/numerics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    /*
     * The log probability that proposeGenotypesFromParentSet() proposes the child's current genotypes: the mean, over
     * its 2T + 1 parent-set options, of the probability of proposing them from that option.
     *
     * An option's probability depends on its observed parents only through their genotypes, and in a clonal cluster
     * most parents share a genome. Parents are grouped by genome (GenomePartition), and each combination of classes,
     * with and without the latent parent, is evaluated once and weighted by the number of options it stands for.
     */
    template<typename InfectionEventImpl, int MaxParentSetSize>
    class ParentSetProposalProbability {
        static_assert(MaxParentSetSize == 1 or MaxParentSetSize == 2, "options are enumerated for at most two observed parents");
        using GeneticsImpl = typename InfectionEventImpl::GeneticsType;

    public:
        // The proposal's rates.
        ParentSetProposalProbability(const double fnRate, const double fpRate) : fnRate_(fnRate), fpRate_(fpRate) {}

        // parentSet is a range of pointers to the observed parents.
        template<typename ParentSet>
        double logProb(const InfectionEventImpl& child, const InfectionEventImpl& latentParent, const ParentSet& parentSet);

    private:
        using Members = std::array<const GeneticsImpl*, MaxParentSetSize + 1>;

        // Log probability of proposing the child's genotypes from the latent parent (if `withLatent`) and the classes'
        // genomes, or -inf.
        double evaluate(bool withLatent, std::span<const unsigned int> classes) const;

        GenomePartition<InfectionEventImpl> partition_{};
        ProposalTarget<GeneticsImpl> child_{};
        std::vector<GeneticsImpl> latentParent_{};
        std::vector<std::vector<GeneticsImpl>> classes_{};// [class][locus]
        std::vector<double> terms_{};

        double fnRate_;
        double fpRate_;
    };

    template<typename InfectionEventImpl, int MaxParentSetSize>
    template<typename ParentSet>
    double ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::logProb(const InfectionEventImpl& child, const InfectionEventImpl& latentParent, const ParentSet& parentSet) {
        partition_.update(parentSet);
        child_.gather(child, fnRate_);
        latentGenotypesAlong(child, latentParent, latentParent_);
        const size_t numClasses = partition_.numClasses();
        classes_.resize(std::max(classes_.size(), numClasses));
        for (std::uint32_t c = 0; c < numClasses; ++c) {
            latentGenotypesAlong(child, partition_.member(c), classes_[c]);
        }

        // The latent parent alone, then every combination of classes with and without it.
        terms_.clear();
        terms_.push_back(evaluate(true, {}));
        double numOptions = 1;
        std::array<unsigned int, 2> combination{};
        auto add = [&](const size_t size, const double count) {
            if (count == 0) {
                return;// a pair within a class of one
            }
            numOptions += 2 * count;
            const double logCount = std::log(count);
            const std::span<const unsigned int> classes(combination.data(), size);
            terms_.push_back(logCount + evaluate(true, classes));
            terms_.push_back(logCount + evaluate(false, classes));
        };
        for (unsigned int c = 0; c < numClasses; ++c) {
            combination[0] = c;
            add(1, partition_.size(c));
        }
        if constexpr (MaxParentSetSize == 2) {
            for (unsigned int c = 0; c < numClasses; ++c) {
                const double size = partition_.size(c);
                combination[0]    = c;
                for (unsigned int d = c; d < numClasses; ++d) {
                    combination[1] = d;
                    add(2, c == d ? size * (size - 1) / 2 : size * partition_.size(d));
                }
            }
        }

        const double maxTerm = *std::max_element(terms_.begin(), terms_.end());
        return utils::logSumExpKnownMax(terms_.begin(), terms_.end(), maxTerm) - std::log(numOptions);
    }

    template<typename InfectionEventImpl, int MaxParentSetSize>
    double ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::evaluate(const bool withLatent, const std::span<const unsigned int> classes) const {
        double logProb = 0;
        Members members{};
        for (size_t l = 0; l < child_.latent.size(); ++l) {
            size_t count = 0;
            if (withLatent) {
                members[count++] = &latentParent_[l];
            }
            for (const auto c : classes) {
                members[count++] = &classes_[c][l];
            }
            const double locusLogProb = logLocusProposalProb(members, count, child_.latent[l], child_.observed[l], child_.presence[l], fpRate_);
            if (locusLogProb == -std::numeric_limits<double>::infinity()) {
                return locusLogProb;
            }
            logProb += locusLogProb;
        }
        return logProb;
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETPROPOSALPROBABILITY_H
