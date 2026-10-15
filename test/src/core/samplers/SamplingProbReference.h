//
// A plain implementation of the proposal probability of RandomAllelesBitSetSampler3
// (see ParentSetGenotypeProposal.h), written from its description rather than for speed: every parent-set option is
// enumerated and each locus' probability is found by enumerating the alleles the chosen parents could contribute.
// Used to check the samplers' calculateSamplingProb(), which reuses work across options and parents. Also the proposal
// density of ParentRelinkSampler, which evaluates each genome once rather than each source.
//

#ifndef TRANSMISSION_NETWORKS_TEST_SAMPLINGPROBREFERENCE_H
#define TRANSMISSION_NETWORKS_TEST_SAMPLINGPROBREFERENCE_H

#include "core/computation/PartialLikelihood.h"
#include "core/datatypes/Genotype.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace transmission_nets::test::reference {
    using transmission_nets::core::datatypes::Genotype;
    using core::computation::Likelihood;

    // The chosen parents' genotypes at a locus: the latent parent if included, then the observed parents at the
    // given indices.
    template<typename InfectionPtr, typename ParentSet, typename Indices, typename LocusPtr>
    std::vector<Genotype> chosenGenotypes(const InfectionPtr& latentParent, const ParentSet& parentSet, const Indices& indices, const bool withLatent, const LocusPtr& locus) {
        std::vector<Genotype> out{};
        if (withLatent) {
            out.push_back(latentParent->latentGenotype(locus)->value());
        }
        for (const auto idx : indices) {
            out.push_back(parentSet.begin()[idx]->latentGenotype(locus)->value());
        }
        return out;
    }

    // Log probability that the proposal draws `child` at one locus from these chosen parents, following its description
    // literally: every combination of one allele picked from each carrying parent, then independent draws of the other
    // carried alleles.
    inline double logLocusProposal(const std::vector<Genotype>& parents, const Genotype& child, const Genotype& observed, const double fnRate, const double fpRate) {
        Genotype carried = child;
        carried.reset();
        std::vector<std::vector<unsigned int>> picks{};
        for (const auto& parent : parents) {
            carried = Genotype::any(carried, parent);
            std::vector<unsigned int> alleles{};
            for (unsigned int i = 0; i < parent.totalAlleles(); ++i) {
                if (parent.allele(i)) {
                    alleles.push_back(i);
                }
            }
            if (not alleles.empty()) {
                picks.push_back(alleles);
            }
        }
        double probability = 0;
        std::vector<size_t> choice(picks.size(), 0);
        while (true) {
            Genotype picked = carried;
            picked.reset();
            double p = 1;
            for (size_t k = 0; k < picks.size(); ++k) {
                picked.set(picks[k][choice[k]], true);
                p /= static_cast<double>(picks[k].size());
            }
            for (unsigned int i = 0; i < carried.totalAlleles(); ++i) {
                const bool present = child.allele(i);
                if (picked.allele(i)) {
                    p *= present ? 1.0 : 0.0;
                } else if (carried.allele(i)) {
                    const double drawn = observed.allele(i) ? 1.0 - fpRate : fnRate;
                    p *= present ? drawn : 1.0 - drawn;
                } else {
                    p *= present ? 0.0 : 1.0;
                }
            }
            probability += p;
            size_t k = 0;
            while (k < picks.size() and ++choice[k] == picks[k].size()) {
                choice[k++] = 0;
            }
            if (k == picks.size()) {
                break;
            }
        }
        return std::log(probability);
    }


    // Log probability of proposing the infection's current genotypes: a uniform mixture over the parent-set options,
    // each the product over loci of logLocusProposal(). A locus with no observation is drawn with every allele
    // unobserved and fnRate 0.5.
    template<int MaxParentSetSize, typename InfectionPtr, typename ParentSet>
    Likelihood proposalProb(const InfectionPtr& infection, const ParentSet& parentSet, const InfectionPtr& latentParent, const double fnRate, const double fpRate) {
        auto optionLogProb = [&](const std::vector<unsigned int>& indices, const bool withLatent) {
            double total = 0;
            for (const auto& locus : infection->loci()) {
                const auto& child = infection->latentGenotype(locus)->value();
                const auto parents = chosenGenotypes(latentParent, parentSet, indices, withLatent, locus);
                if (infection->observedGenotype().contains(locus)) {
                    total += logLocusProposal(parents, child, infection->observedGenotype(locus)->value(), fnRate, fpRate);
                } else {
                    Genotype nothingObserved = child;
                    nothingObserved.reset();
                    total += logLocusProposal(parents, child, nothingObserved, 0.5, fpRate);
                }
            }
            return total;
        };

        std::vector<Likelihood> options{};
        for (size_t k = 1; k <= MaxParentSetSize and k <= parentSet.size(); ++k) {
            core::utils::generators::CombinationIndicesGenerator combinations(parentSet.size(), k);
            while (not combinations.completed) {
                const std::vector<unsigned int> indices(combinations.curr.begin(), combinations.curr.end());
                options.push_back(optionLogProb(indices, true));
                options.push_back(optionLogProb(indices, false));
                combinations.next();
            }
        }
        options.push_back(optionLogProb({}, true));

        const Likelihood maxLogProb = *std::max_element(options.begin(), options.end());
        return core::utils::logSumExpKnownMax(options.begin(), options.end(), maxLogProb) - std::log(static_cast<double>(options.size()));
    }

    // Log density of ParentRelinkSampler's proposal at the infection's current duration and genotypes: a uniform
    // mixture over the sources (allowed parents, then the latent parent) whose window of durations holds positive
    // shape mass, each the shape restricted to that window times the probability of drawing every locus from the
    // source alone. A duration d falls in day round(d), and a source's window ends where the infection would no longer
    // be infected after it.
    template<typename InfectionPtr>
    double relinkProposalDensity(const InfectionPtr& infection, const InfectionPtr& latentParent, const std::vector<InfectionPtr>& allowedParents,
                                 const double minDuration, const double maxDuration, const std::vector<double>& shape, const double fnRate, const double fpRate) {
        auto mass = [&](const double first, const double second) {
            double total = 0;
            for (size_t day = 0; day < shape.size(); ++day) {
                const double overlap = std::min(second, day + 0.5) - std::max(first, day - 0.5);
                if (overlap > 0) {
                    total += shape[day] * overlap;
                }
            }
            return total;
        };
        auto genomeLogProb = [&](const InfectionPtr& source) {
            double total = 0;
            for (const auto& locus : infection->loci()) {
                const auto& child = infection->latentGenotype(locus)->value();
                const std::vector<Genotype> parents{source->latentGenotype(locus)->value()};
                if (infection->observedGenotype().contains(locus)) {
                    total += logLocusProposal(parents, child, infection->observedGenotype(locus)->value(), fnRate, fpRate);
                } else {
                    Genotype nothingObserved = child;
                    nothingObserved.reset();
                    total += logLocusProposal(parents, child, nothingObserved, 0.5, fpRate);
                }
            }
            return total;
        };

        const double duration = infection->infectionDuration()->value();
        const long day        = std::lround(duration);
        if (day < 0 or day >= static_cast<long>(shape.size()) or shape[day] <= 0) {
            return -std::numeric_limits<double>::infinity();
        }
        std::vector<double> terms{};
        size_t feasible = 0;
        for (size_t s = 0; s <= allowedParents.size(); ++s) {
            const bool latent   = s == allowedParents.size();
            const double first  = minDuration;
            const double second = latent ? maxDuration : std::min(maxDuration, infection->observationTime()->value() - allowedParents[s]->infectionTime());
            if (second <= first) {
                continue;
            }
            const double sourceMass = mass(first, second);
            if (sourceMass <= 0) {
                continue;
            }
            ++feasible;
            if (duration >= first and duration < second) {
                terms.push_back(std::log(shape[day] / sourceMass) + genomeLogProb(latent ? latentParent : allowedParents[s]));
            }
        }
        if (terms.empty()) {
            return -std::numeric_limits<double>::infinity();
        }
        const double maxTerm = *std::max_element(terms.begin(), terms.end());
        return core::utils::logSumExpKnownMax(terms.begin(), terms.end(), maxTerm) - std::log(static_cast<double>(feasible));
    }

}// namespace transmission_nets::test::reference

#endif//TRANSMISSION_NETWORKS_TEST_SAMPLINGPROBREFERENCE_H
