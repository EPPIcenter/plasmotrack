//
// Checks that RandomAllelesBitSetSampler3 and ParentRelinkSampler compute the same proposal
// probabilities as the plain implementations in SamplingProbReference.h.
//

#include "gtest/gtest.h"

#include "SamplingProbReference.h"

#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/specialized/ParentRelinkSampler.h"

#include <boost/random.hpp>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::samplers;

namespace {
    // Same configuration as the model (src/impl/model/Model/config.h).
    constexpr int MAX_ALLELES = 128;
    constexpr int MAX_PARENTS = 2;
    constexpr int MAX_COI     = 20;

    using GeneticsImpl   = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = Infection<GeneticsImpl>;
    using Engine         = boost::random::mt19937;

    // calculateSamplingProb only reads the Candidate parents; the target is never used (update() is not called).
    struct FixedOrdering {
        ParentSet<InfectionEvent> parents{};
        [[nodiscard]] const ParentSet<InfectionEvent>& candidateParents(const std::shared_ptr<InfectionEvent>&) const { return parents; }
    };
    struct UnusedTarget {
        double value() { return 0; }
        [[nodiscard]] bool isDirty() const { return false; }
    };

    using Sampler3 = genetics::RandomAllelesBitSetSampler3<UnusedTarget, Engine, InfectionEvent, GeneticsImpl, FixedOrdering, MAX_PARENTS, MAX_COI>;
    using Relink   = specialized::ParentRelinkSampler<UnusedTarget, Engine, InfectionEvent>;

    std::string randomGenotype(Engine& rng, const unsigned int numAlleles) {
        boost::random::uniform_01<> unif;
        boost::random::uniform_int_distribution<> pick(0, static_cast<int>(numAlleles) - 1);
        std::string bits(numAlleles, '0');
        bits[pick(rng)] = '1';
        if (unif(rng) < 0.3) {
            bits[pick(rng)] = '1';
        }
        return bits;
    }
}// namespace

TEST(SamplingProbTest, MatchesReferenceImplementationsExactly) {
    Engine rng(20261005);
    boost::random::uniform_01<> unif;

    const std::vector<unsigned int> alleleCounts{2, 3, 4, 5, 8, 13};
    std::vector<std::shared_ptr<Locus>> loci;
    for (size_t l = 0; l < alleleCounts.size(); ++l) {
        loci.push_back(std::make_shared<Locus>("L" + std::to_string(l), alleleCounts[l]));
    }

    // Clonal founders with occasional changes, and some loci left unobserved (latent only).
    constexpr int numFounders   = 3;
    constexpr int numInfections = 16;
    std::vector<std::vector<std::string>> founders(numFounders);
    for (auto& founder : founders) {
        for (const auto count : alleleCounts) {
            founder.push_back(randomGenotype(rng, count));
        }
    }
    std::vector<std::shared_ptr<InfectionEvent>> infections;
    std::vector<std::shared_ptr<InfectionEvent>> latentParents;
    for (int i = 0; i < numInfections; ++i) {
        auto inf = std::make_shared<InfectionEvent>("inf" + std::to_string(i), 100.0 + i, true);
        for (size_t l = 0; l < loci.size(); ++l) {
            const auto observed = unif(rng) < 0.15 ? randomGenotype(rng, alleleCounts[l]) : founders[i % numFounders][l];
            const auto latent   = unif(rng) < 0.1 ? randomGenotype(rng, alleleCounts[l]) : observed;
            if (unif(rng) < 0.1) {
                inf->addLatentGenetics(loci[l], latent);
            } else {
                inf->addGenetics(loci[l], observed, latent);
            }
        }
        infections.push_back(inf);
        latentParents.push_back(std::make_shared<InfectionEvent>(*inf));
    }

    // The reference and the samplers sum the same terms in different orders, so finite values agree to rounding.
    auto expectSame = [](const double expected, const double actual, const std::string& label) {
        if (std::isfinite(expected)) {
            EXPECT_NEAR(expected, actual, 1e-9 * std::max(1.0, std::abs(expected))) << label;
        } else {
            EXPECT_EQ(expected, actual) << label;
        }
    };

    // Durations spread over the duration shape's range, so that the relink sampler's sources have different windows.
    std::vector<double> shape(60);
    for (size_t day = 1; day < shape.size(); ++day) {
        shape[day] = 1.0 + static_cast<double>(day);
    }
    for (const auto& infection : infections) {
        infection->infectionDuration()->saveState(1);
        infection->infectionDuration()->setValue(1.0 + 58.0 * unif(rng));
        infection->infectionDuration()->acceptState();
    }

    // Samplers are kept across rounds, as in a run: each keeps its partition of parents by genome between
    // evaluations, and must notice when their genomes change.
    std::vector<std::shared_ptr<FixedOrdering>> parentSets;
    std::vector<std::unique_ptr<Sampler3>> samplers3;
    std::vector<std::unique_ptr<Relink>> relinks;
    for (int c = 0; c < numInfections; ++c) {
        std::vector<std::shared_ptr<InfectionEvent>> allowed;
        for (int p = 0; p < numInfections; ++p) {
            if (p != c) {
                allowed.push_back(infections[p]);
            }
        }
        parentSets.push_back(std::make_shared<FixedOrdering>());
        samplers3.push_back(std::make_unique<Sampler3>(infections[c], parentSets[c], latentParents[c], nullptr, std::make_shared<Engine>(1), 0.1, 0.1));
        relinks.push_back(std::make_unique<Relink>(infections[c], latentParents[c], allowed, nullptr, std::make_shared<Engine>(1), 1.0, 59.0, shape, 0.05, 0.05));
    }

    int comparisons = 0;
    int finite      = 0;
    int finiteRelink = 0;
    auto check = [&](const int c, const std::string& where) {
        std::vector<std::shared_ptr<InfectionEvent>> allowed;
        for (int p = 0; p < numInfections; ++p) {
            if (p != c) {
                allowed.push_back(infections[p]);
            }
        }
        const double expected3 = transmission_nets::test::reference::proposalProb<GeneticsImpl, MAX_PARENTS>(infections[c], parentSets[c]->parents, latentParents[c], 0.1, 0.1);
        const double expectedR = transmission_nets::test::reference::relinkProposalDensity<GeneticsImpl>(infections[c], latentParents[c], allowed, 1.0, 59.0, shape, 0.05, 0.05);
        const auto label = where + ", child " + std::to_string(c) + ", parents " + std::to_string(parentSets[c]->parents.size());
        expectSame(expected3, samplers3[c]->calculateSamplingProb(), "Sampler3 " + label);
        expectSame(expectedR, relinks[c]->logProposalDensity(), "Relink " + label);
        comparisons += 2;
        finite += std::isfinite(expected3) + std::isfinite(expectedR);
        finiteRelink += std::isfinite(expectedR);
    };

    for (int round = 0; round < 3; ++round) {
        for (int c = 0; c < numInfections; ++c) {
            const auto where = "round " + std::to_string(round);
            // The previous round's parent set, whose members' genomes have since changed.
            if (round > 0) {
                check(c, where + ", same parents");
            }
            // Parent sets of every size from none to all other infections.
            parentSets[c]->parents.clear();
            for (int p = 0; p < numInfections; ++p) {
                if (p != c and unif(rng) < (round + 1) / 3.0) {
                    parentSets[c]->parents.insert(infections[p]);
                }
            }
            check(c, where);

            // Samplers evaluate before and after each proposal, so the same sampler must not reuse earlier results.
            const auto& locus = loci[(c + round) % loci.size()];
            auto genotype     = infections[c]->latentGenotype(locus);
            genotype->saveState(1);
            genotype->setValue(GeneticsImpl(randomGenotype(rng, locus->totalAlleles())));
            genotype->acceptState();
            check(c, where + ", re-evaluated");
        }

        // Change some latent genotypes between rounds, and give some infections another's whole genome so that
        // genome classes have several members.
        for (int i = 0; i < numInfections; ++i) {
            const auto& target = i % 2 == 0 ? infections[i] : latentParents[i];
            const auto& locus  = loci[(i + round) % loci.size()];
            auto genotype      = target->latentGenotype(locus);
            genotype->saveState(1);
            genotype->setValue(GeneticsImpl(randomGenotype(rng, locus->totalAlleles())));
            genotype->acceptState();
        }
        for (int i = round; i + 1 < numInfections; i += 3) {
            for (const auto& locus : loci) {
                auto genotype = infections[i + 1]->latentGenotype(locus);
                genotype->saveState(1);
                genotype->setValue(infections[i]->latentGenotype(locus)->value());
                genotype->acceptState();
            }
        }
    }

    EXPECT_EQ(comparisons, 2 * numInfections * 8);
    // Guard against a test that only ever compares -inf. Many of these random parent sets are impossible: every
    // chosen parent must share an allele with the child at every locus, observed or not, and genotype changes
    // accumulate over the rounds while the latent parents keep their starting genotypes.
    EXPECT_GT(finite, comparisons / 4);
    EXPECT_GT(finiteRelink, comparisons / 8);
}
