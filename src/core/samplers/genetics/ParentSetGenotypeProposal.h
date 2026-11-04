//
// The genotype proposal of RandomAllelesBitSetSampler3.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETGENOTYPEPROPOSAL_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETGENOTYPEPROPOSAL_H

#include "core/containers/UidLess.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/random.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    // At a locus with no observation, every allele carried by a chosen parent (beyond the one it contributes) is
    // present with this probability. Such loci are drawn as if every allele were unobserved, with fnRate set to it.
    inline constexpr double UNOBSERVED_LOCUS_PRESENCE = 0.5;

    /*
     * The proposal at one locus, given the chosen parents' genotypes there. Each chosen parent carrying alleles
     * contributes one of them, picked uniformly, so that it could have transmitted to the child. Every other allele
     * carried by a chosen parent is then present with probability 1 - fpRate if observed in the child and fnRate if
     * not, and every allele no chosen parent carries is absent.
     */
    template<typename GeneticsImpl, size_t N, typename Engine>
    GeneticsImpl proposeLocus(const std::array<const GeneticsImpl*, N>& parents, const size_t numParents, const GeneticsImpl& observed, Engine& rng, const double fnRate, const double fpRate) {
        auto carried = observed;
        carried.reset();
        auto proposal = carried;
        for (size_t p = 0; p < numParents; ++p) {
            carried = GeneticsImpl::any(carried, *parents[p]);
            if (const unsigned int count = parents[p]->totalPositiveCount(); count > 0) {
                boost::random::uniform_int_distribution<unsigned int> pick(0, count - 1);
                unsigned int skip = pick(rng);
                for (unsigned int i = 0; i < parents[p]->totalAlleles(); ++i) {
                    if (parents[p]->allele(i) and skip-- == 0) {
                        proposal.set(i, true);
                        break;
                    }
                }
            }
        }
        boost::random::uniform_01<> unif;
        for (unsigned int i = 0; i < carried.totalAlleles(); ++i) {
            if (carried.allele(i) and not proposal.allele(i)) {
                proposal.set(i, unif(rng) < (observed.allele(i) ? 1.0 - fpRate : fnRate));
            }
        }
        return proposal;
    }

    /*
     * Log probability that proposeLocus() gives `child`, or -inf. A sum over the alleles the carrying parents could
     * have contributed (each one in the child): every carried allele is drawn independently, except that a
     * contributed allele is present with probability 1 instead of its draw's.
     */
    template<typename GeneticsImpl, size_t N>
    double logLocusProposalProb(const std::array<const GeneticsImpl*, N>& parents, const size_t numParents, const GeneticsImpl& child, const GeneticsImpl& observed, const double fnRate, const double fpRate) {
        auto carried = child;
        carried.reset();
        for (size_t p = 0; p < numParents; ++p) {
            carried = GeneticsImpl::any(carried, *parents[p]);
        }
        if (GeneticsImpl::falsePositiveCount(carried, child) > 0) {
            return -std::numeric_limits<double>::infinity();// the child has an allele no chosen parent carries
        }

        // Every carried allele drawn independently.
        const auto presentObserved    = GeneticsImpl::truePositiveCount(child, observed);
        const auto presentUnobserved  = GeneticsImpl::falseNegativeCount(child, observed);
        const auto absentObserved     = GeneticsImpl::truePositiveCount(carried, observed) - presentObserved;
        const auto absentUnobserved   = GeneticsImpl::falseNegativeCount(carried, observed) - presentUnobserved;
        const double logIndependent   = presentObserved * std::log(1.0 - fpRate) + absentObserved * std::log(fpRate) +
                                        presentUnobserved * std::log(fnRate) + absentUnobserved * std::log(1.0 - fnRate);

        // Each carrying parent's possible contributions: its alleles that are in the child.
        auto nextAllele = [](const GeneticsImpl& alleles, unsigned int from) {
            while (from < alleles.totalAlleles() and not alleles.allele(from)) {
                ++from;
            }
            return from;
        };
        std::array<const GeneticsImpl*, N> carrying{};
        std::array<GeneticsImpl, N> inChild{};
        std::array<unsigned int, N> choice{};// the contributed allele, as an index
        size_t numCarrying = 0;
        for (size_t p = 0; p < numParents; ++p) {
            if (parents[p]->totalPositiveCount() == 0) {
                continue;
            }
            inChild[numCarrying] = GeneticsImpl::shared(*parents[p], child);
            if (inChild[numCarrying].totalPositiveCount() == 0) {
                return -std::numeric_limits<double>::infinity();// this parent contributed nothing
            }
            choice[numCarrying]   = nextAllele(inChild[numCarrying], 0);
            carrying[numCarrying] = parents[p];
            ++numCarrying;
        }

        // Sum over every combination of contributions of the probability of picking it, times 1 / (draw probability)
        // for each distinct contributed allele, whose draw it replaces.
        double total = 0;
        while (true) {
            double weight = 1;
            std::array<unsigned int, N> contributed{};
            size_t numContributed = 0;
            for (size_t p = 0; p < numCarrying; ++p) {
                weight /= carrying[p]->totalPositiveCount();
                const unsigned int allele = choice[p];
                if (std::find(contributed.begin(), contributed.begin() + numContributed, allele) == contributed.begin() + numContributed) {
                    contributed[numContributed++] = allele;
                    weight /= observed.allele(allele) ? 1.0 - fpRate : fnRate;
                }
            }
            total += weight;

            size_t p = 0;
            for (; p < numCarrying; ++p) {
                choice[p] = nextAllele(inChild[p], choice[p] + 1);
                if (choice[p] < inChild[p].totalAlleles()) {
                    break;
                }
                choice[p] = nextAllele(inChild[p], 0);
            }
            if (p == numCarrying) {
                break;
            }
        }
        return logIndependent + std::log(total);
    }

    /*
     * An infection's genotypes at each of its loci as the proposal sees them, in latentGenotype() order: the latent
     * genotype, the observed one (empty where nothing was observed), and the presence rate proposeLocus() uses for
     * carried alleles the observation lacks (fnRate, or UNOBSERVED_LOCUS_PRESENCE where nothing was observed).
     */
    template<typename GeneticsImpl>
    struct ProposalTarget {
        std::vector<GeneticsImpl> latent{};
        std::vector<GeneticsImpl> observed{};
        std::vector<double> presence{};

        template<typename InfectionEventImpl>
        void gather(const InfectionEventImpl& infection, const double fnRate) {
            latent.clear();
            observed.clear();
            presence.clear();
            // Both maps are ordered by locus, so the observations are merged in rather than looked up.
            const auto& observations = infection.observedGenotype();
            auto observation         = observations.begin();
            for (const auto& [locus, genotype] : infection.latentGenotype()) {
                latent.push_back(genotype->value());
                while (observation != observations.end() and containers::UidLess{}(observation->first, locus)) {
                    ++observation;
                }
                if (observation != observations.end() and observation->first == locus) {
                    observed.push_back(observation->second->value());
                    presence.push_back(fnRate);
                } else {
                    observed.push_back(genotype->value());
                    observed.back().reset();
                    presence.push_back(UNOBSERVED_LOCUS_PRESENCE);
                }
            }
        }
    };

    // `other`'s latent genotypes at `infection`'s loci, in infection.latentGenotype() order.
    template<typename InfectionEventImpl, typename GeneticsImpl>
    void latentGenotypesAlong(const InfectionEventImpl& infection, const InfectionEventImpl& other, std::vector<GeneticsImpl>& out) {
        out.clear();
        const auto& theirs = other.latentGenotype();
        auto it            = theirs.begin();
        for (const auto& [locus, genotype] : infection.latentGenotype()) {
            // Infections normally have the same loci, kept in the same order, so this lookup is rarely needed.
            if (it == theirs.end() or it->first != locus) {
                it = theirs.find(locus);
                if (it == theirs.end()) {
                    throw std::out_of_range("infection lacks a locus of the infection proposed for");
                }
            }
            out.push_back(it->second->value());
            ++it;
        }
    }

    /*
     * `latentParent` here is the infection's latent-parent placeholder, whose genotypes are a fixed anchor (the observed
     * genotypes, or every allele where nothing was observed; see State). The latent parent itself is summed out of the
     * model, so the anchor only shapes the proposal, which stays valid because the anchor never changes.
     *
     * Proposes an infection's latent genotypes from one parent-set option, chosen uniformly among 2T + 1, T being the
     * number of combinations of up to MaxParentSetSize observed parents: each combination with or without the latent
     * parent, or the latent parent alone. Each locus is then drawn by proposeLocus() from the chosen parents, a locus
     * with no observation as described at UNOBSERVED_LOCUS_PRESENCE. Missing loci matter: a child can only take an
     * observed parent if every locus is compatible with it.
     *
     * The samplers' calculateSamplingProb() is the probability of this proposal, so the two must change together.
     * Every locus is saved under savedStateId, changed or not.
     */
    template<int MaxParentSetSize, typename InfectionPtr, typename ParentSet, typename Engine>
    void proposeGenotypesFromParentSet(const InfectionPtr& infection, const InfectionPtr& latentParent, const ParentSet& parentSet, Engine& rng, const double fnRate, const double fpRate, const int savedStateId) {
        const size_t numParents = parentSet.size();
        size_t numCombinations  = 0;
        for (size_t k = 1; k <= MaxParentSetSize and k <= numParents; ++k) {
            numCombinations += static_cast<size_t>(boost::math::binomial_coefficient<double>(numParents, k));
        }

        boost::random::uniform_int_distribution<size_t> pickOption(0, 2 * numCombinations);
        const size_t option     = pickOption(rng);
        const bool latentOnly   = option == 2 * numCombinations;
        const bool withLatent   = latentOnly or option % 2 == 0;
        std::vector<unsigned int> chosen{};
        if (not latentOnly) {
            size_t combination = option / 2;
            utils::generators::CombinationIndicesGenerator combinations;
            for (size_t k = 1; k <= MaxParentSetSize and k <= numParents; ++k) {
                const auto count = static_cast<size_t>(boost::math::binomial_coefficient<double>(numParents, k));
                if (combination < count) {
                    combinations.reset(numParents, k);
                    for (size_t skipped = 0; skipped < combination; ++skipped) {
                        combinations.next();
                    }
                    chosen.assign(combinations.curr.begin(), combinations.curr.end());
                    break;
                }
                combination -= count;
            }
        }

        for (const auto& locus : infection->loci()) {
            auto genotype = infection->latentGenotype(locus);
            genotype->saveState(savedStateId);

            using Genetics = std::remove_cvref_t<decltype(genotype->value())>;
            std::array<const Genetics*, MaxParentSetSize + 1> parents{};
            size_t numChosen = 0;
            if (withLatent) {
                parents[numChosen++] = &latentParent->latentGenotype(locus)->value();
            }
            for (const auto idx : chosen) {
                parents[numChosen++] = &parentSet.begin()[idx]->latentGenotype(locus)->value();
            }
            Genetics proposal;
            if (infection->observedGenotype().contains(locus)) {
                proposal = proposeLocus(parents, numChosen, infection->observedGenotype(locus)->value(), rng, fnRate, fpRate);
            } else {
                auto nothingObserved = genotype->value();
                nothingObserved.reset();
                proposal = proposeLocus(parents, numChosen, nothingObserved, rng, UNOBSERVED_LOCUS_PRESENCE, fpRate);
            }
            if (proposal != genotype->value()) {
                genotype->setValue(proposal);
            }
        }
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETGENOTYPEPROPOSAL_H
