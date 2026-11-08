//
// The proposal probability of RandomAllelesBitSetSampler3.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETPROPOSALPROBABILITY_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETPROPOSALPROBABILITY_H

#include "core/datatypes/Genotype.h"
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
     * its 2T + 1 parent-set options, of the probability of proposing them from that option. With option weights, the
     * weighted mean instead.
     *
     * An option's probability depends on its observed parents only through their genotypes, and in a clonal cluster
     * most parents share a genome. Parents are grouped by genome (GenomePartition), and each combination of classes,
     * with and without the latent parent, is evaluated once and weighted by the options it stands for.
     */
    template<typename InfectionEventImpl, int MaxParentSetSize>
    class ParentSetProposalProbability {
        static_assert(MaxParentSetSize == 1 or MaxParentSetSize == 2, "options are enumerated for at most two observed parents");
        using Genotype = core::datatypes::Genotype;

    public:
        // Per group of options, the log probability of proposing the child's genotypes from any one of them.
        using Evaluation = std::vector<double>;

        // The proposal's rates.
        ParentSetProposalProbability(const double fnRate, const double fpRate) : fnRate_(fnRate), fpRate_(fpRate) {}

        void setRates(const double fnRate, const double fpRate) {
            fnRate_ = fnRate;
            fpRate_ = fpRate;
        }

        // parentSet is a range of pointers to the observed parents.
        template<typename ParentSet>
        double logProb(const InfectionEventImpl& child, const InfectionEventImpl& latentParent, const ParentSet& parentSet) {
            evaluate(child, latentParent, parentSet, scratch_);
            return combine(scratch_, {});
        }

        // Evaluates the child's current genotypes against every group of options. The parent set's genomes must not
        // change between this and combine().
        // Only the loci `include` accepts are evaluated, as when only they are proposed.
        template<typename ParentSet, typename Include = AllLoci>
        void evaluate(const InfectionEventImpl& child, const InfectionEventImpl& latentParent, const ParentSet& parentSet, Evaluation& out, const Include& include = {});

        // Re-evaluates the child's genotypes against the groups already built by the preceding evaluate() call when
        // only the child's latent genotypes have changed.
        void evaluateChildOnly(const InfectionEventImpl& child, Evaluation& out);

        // The log proposal probability of an evaluated state. `optionProbs`, if not empty, gives each option's
        // probability in parent-set-distribution order (see parentSetIndexOfOption()); otherwise options are uniform.
        double combine(const Evaluation& evaluation, std::span<const double> optionProbs);

    private:
        using Members = std::array<const core::datatypes::Genotype*, MaxParentSetSize + 1>;
        static constexpr size_t NO_GROUP = std::numeric_limits<size_t>::max();

        struct Group {
            bool withLatent;
            std::array<unsigned int, 2> classes;
            size_t size;
            double count;// options it stands for
        };

        // Log probability of proposing the child's genotypes from the latent parent (if `withLatent`) and the classes'
        // genomes, or -inf.
        double evaluateGroup(bool withLatent, std::span<const unsigned int> classes) const;

        GenomePartition<InfectionEventImpl> partition_{};
        ProposalTarget child_{};
        std::vector<core::datatypes::Genotype> latentParent_{};
        std::vector<std::vector<core::datatypes::Genotype>> classes_{};// [class][locus]
        std::vector<size_t> includedLoci_{};

        std::vector<Group> groups_{};
        // The group of the options of class c alone with the latent parent (the one without it follows), and of the
        // pair of classes c <= d at [c * numClasses + d].
        std::vector<size_t> singleGroup_{};
        std::vector<size_t> pairGroup_{};
        size_t numParents_ = 0;
        double numOptions_ = 1;

        Evaluation scratch_{};
        std::vector<double> groupWeights_{};
        std::vector<double> terms_{};

        double fnRate_;
        double fpRate_;
    };

    template<typename InfectionEventImpl, int MaxParentSetSize>
    template<typename ParentSet, typename Include>
    void ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::evaluate(const InfectionEventImpl& child, const InfectionEventImpl& latentParent, const ParentSet& parentSet, Evaluation& out, const Include& include) {
        partition_.update(parentSet);
        numParents_ = parentSet.size();
        child_.gather(child, fnRate_);
        includedLoci_.clear();
        size_t l = 0;
        for (const auto& [locus, genotype] : child.latentGenotype()) {
            if (include(locus)) {
                includedLoci_.push_back(l);
            }
            ++l;
        }
        latentGenotypesAlong(child, latentParent, latentParent_);
        const size_t numClasses = partition_.numClasses();
        classes_.resize(std::max(classes_.size(), numClasses));
        for (std::uint32_t c = 0; c < numClasses; ++c) {
            latentGenotypesAlong(child, partition_.member(c), classes_[c]);
        }

        // The latent parent alone, then every combination of classes with and without it.
        groups_.clear();
        singleGroup_.assign(numClasses, NO_GROUP);
        pairGroup_.assign(numClasses * numClasses, NO_GROUP);
        groups_.push_back({.withLatent = true, .classes = {}, .size = 0, .count = 1});
        numOptions_ = 1;
        auto add = [&](const std::array<unsigned int, 2> classes, const size_t size, const double count) {
            if (count == 0) {
                return NO_GROUP;// a pair within a class of one
            }
            numOptions_ += 2 * count;
            groups_.push_back({.withLatent = true, .classes = classes, .size = size, .count = count});
            groups_.push_back({.withLatent = false, .classes = classes, .size = size, .count = count});
            return groups_.size() - 2;
        };
        for (unsigned int c = 0; c < numClasses; ++c) {
            singleGroup_[c] = add({c, 0}, 1, partition_.size(c));
        }
        if constexpr (MaxParentSetSize == 2) {
            for (unsigned int c = 0; c < numClasses; ++c) {
                const double size = partition_.size(c);
                for (unsigned int d = c; d < numClasses; ++d) {
                    pairGroup_[c * numClasses + d] = add({c, d}, 2, c == d ? size * (size - 1) / 2 : size * partition_.size(d));
                }
            }
        }

        out.clear();
        for (const auto& group : groups_) {
            out.push_back(evaluateGroup(group.withLatent, std::span<const unsigned int>(group.classes.data(), group.size)));
        }
    }

    template<typename InfectionEventImpl, int MaxParentSetSize>
    void ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::evaluateChildOnly(const InfectionEventImpl& child, Evaluation& out) {
        child_.gather(child, fnRate_);
        out.clear();
        for (const auto& group : groups_) {
            out.push_back(evaluateGroup(group.withLatent, std::span<const unsigned int>(group.classes.data(), group.size)));
        }
    }

    template<typename InfectionEventImpl, int MaxParentSetSize>
    double ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::combine(const Evaluation& evaluation, const std::span<const double> optionProbs) {
        terms_.clear();
        if (optionProbs.empty()) {
            for (size_t g = 0; g < groups_.size(); ++g) {
                terms_.push_back(std::log(groups_[g].count) + evaluation[g]);
            }
            const double maxTerm = *std::max_element(terms_.begin(), terms_.end());
            return utils::logSumExpKnownMax(terms_.begin(), terms_.end(), maxTerm) - std::log(numOptions_);
        }

        // Each group's weight, summed over its options in parent-set-distribution order: the latent parent alone,
        // then each single parent and each pair, without and with the latent parent.
        groupWeights_.assign(groups_.size(), 0.0);
        groupWeights_[0] = optionProbs[0];
        const size_t numClasses = partition_.numClasses();
        size_t index            = 1;
        auto add                = [&](const size_t group) {
            groupWeights_[group + 1] += optionProbs[index];
            groupWeights_[group] += optionProbs[index + 1];
            index += 2;
        };
        for (size_t i = 0; i < numParents_; ++i) {
            add(singleGroup_[partition_.classOf(i)]);
        }
        if constexpr (MaxParentSetSize == 2) {
            for (size_t i = 0; i + 1 < numParents_; ++i) {
                for (size_t j = i + 1; j < numParents_; ++j) {
                    const auto c = std::min(partition_.classOf(i), partition_.classOf(j));
                    const auto d = std::max(partition_.classOf(i), partition_.classOf(j));
                    add(pairGroup_[c * numClasses + d]);
                }
            }
        }
        for (size_t g = 0; g < groups_.size(); ++g) {
            terms_.push_back(std::log(groupWeights_[g]) + evaluation[g]);
        }
        const double maxTerm = *std::max_element(terms_.begin(), terms_.end());
        return utils::logSumExpKnownMax(terms_.begin(), terms_.end(), maxTerm);
    }

    template<typename InfectionEventImpl, int MaxParentSetSize>
    double ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize>::evaluateGroup(const bool withLatent, const std::span<const unsigned int> classes) const {
        double logProb = 0;
        Members members{};
        for (const size_t l : includedLoci_) {
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
