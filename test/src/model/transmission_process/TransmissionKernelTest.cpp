//
// The Transmission kernel from plain genotypes: the locus term on cases small enough to work by hand, the strain-count
// table and the sum over strain counts, and the latent locus table against the plain computation in
// MarginalLatentParentReference.h.
//

#include "gtest/gtest.h"

#include "MarginalLatentParentReference.h"

#include "core/datatypes/Alleles.h"
#include "model/transmission_process/TransmissionKernel.h"

#include <boost/random.hpp>

#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <vector>

using namespace transmission_nets::model::transmission_process;
using transmission_nets::core::datatypes::AllelesBitSet;

namespace {
    constexpr int MAX_PARENTS = 2;
    constexpr int MAX_STRAINS = 12;
    constexpr int MAX_COI     = 20;
    constexpr double NEG_INF  = -std::numeric_limits<double>::infinity();

    using Genetics = AllelesBitSet<128>;
    using Kernel   = TransmissionKernel<MAX_PARENTS, MAX_STRAINS, MAX_COI, Genetics>;
    using Variant  = Kernel::Variant;
    using Engine   = boost::random::mt19937;

    std::vector<const Genetics*> pointers(const std::vector<Genetics>& genotypes) {
        std::vector<const Genetics*> out;
        for (const auto& g : genotypes) {
            out.push_back(&g);
        }
        return out;
    }

    Genetics fromMask(const unsigned int numAlleles, const unsigned int mask) {
        Genetics g(static_cast<int>(numAlleles));
        g.reset();
        for (unsigned int i = 0; i < numAlleles; ++i) {
            if (mask >> i & 1u) {
                g.set(i);
            }
        }
        return g;
    }
    unsigned int randomMask(Engine& rng, const unsigned int numAlleles, const double extra) {
        boost::random::uniform_int_distribution<unsigned int> pick(0, numAlleles - 1);
        boost::random::uniform_01<> u;
        unsigned int mask = 1u << pick(rng);
        while (u(rng) < extra) {
            mask |= 1u << pick(rng);
        }
        return mask;
    }
    std::vector<double> randomFrequencies(Engine& rng, const unsigned int numAlleles) {
        boost::random::uniform_01<> u;
        std::vector<double> f(numAlleles);
        double s = 0;
        for (auto& x : f) {
            x = 0.05 + u(rng);
            s += x;
        }
        for (auto& x : f) {
            x /= s;
        }
        return f;
    }
}// namespace

TEST(TransmissionKernelTest, AParentCarryingTheChildsOnlyAlleleTransmitsItWithCertainty) {
    Kernel kernel;
    const Genetics child("10"), parent("10");
    const std::vector<const Genetics*> parents{&parent};
    const auto term = kernel.locusTerm(Variant::ObservedParents, child, parents, nullptr);
    ASSERT_FALSE(term.impossible);
    for (const double t : term.terms) {
        EXPECT_DOUBLE_EQ(t, 0.0);
    }
}

TEST(TransmissionKernelTest, TwoChildAllelesNeedTwoStrainsAndEachStrainPicksOneAtRandom) {
    Kernel kernel;
    const Genetics child("11"), parent("11");
    const std::vector<const Genetics*> parents{&parent};
    const auto term = kernel.locusTerm(Variant::ObservedParents, child, parents, nullptr);
    ASSERT_FALSE(term.impossible);
    EXPECT_EQ(term.terms[0], NEG_INF);           // one strain cannot carry both
    EXPECT_NEAR(term.terms[1], std::log(0.5), 1e-12);// two strains: both alleles present unless both draws agree
}

TEST(TransmissionKernelTest, AParentLackingAChildAlleleIsImpossible) {
    Kernel kernel;
    const Genetics child("11"), parent("10"), latent("01"), unrelated("001");
    const std::vector<const Genetics*> parents{&parent};
    EXPECT_TRUE(kernel.locusTerm(Variant::ObservedParents, child, parents, nullptr).impossible);
    // A latent parent sharing no allele with the child cannot be one of its parents.
    EXPECT_TRUE(kernel.locusTerm(Variant::ObservedAndLatentParents, child, parents, &unrelated).impossible);
    // One that supplies the missing allele makes the set possible.
    EXPECT_FALSE(kernel.locusTerm(Variant::ObservedAndLatentParents, child, parents, &latent).impossible);
    EXPECT_FALSE(kernel.locusTerm(Variant::LatentParentOnly, child, {}, &child).impossible);
}

TEST(TransmissionKernelTest, StrainCountTableRowsAreDistributionsOverCountsAtLeastTheParentCount) {
    const auto table = Kernel::strainCountTable(1.7);
    for (int n = 1; n <= MAX_PARENTS + 1; ++n) {
        double total = 0;
        for (int m = 1; m <= MAX_STRAINS; ++m) {
            const double logP = Kernel::probNumStrains(table, m, n);
            if (m < n) {
                EXPECT_EQ(logP, NEG_INF) << "n=" << n << " m=" << m;
            } else {
                EXPECT_GT(logP, NEG_INF) << "n=" << n << " m=" << m;
                total += std::exp(logP);
            }
        }
        EXPECT_NEAR(total, 1.0, 1e-12) << "n=" << n;
    }
}

TEST(TransmissionKernelTest, StrainSumsExcludeCountsBelowTheNumberOfParents) {
    // Every strain carries the child's only allele, so each locus term is log(1) = 0 and the sum over strain counts
    // reduces to the strain-count distribution over the possible counts, which is 1. Counting the impossible
    // single-strain case for two parents used to add exactly log(2).
    const auto table = Kernel::strainCountTable(2.0);
    const Kernel::StrainSums sums{};
    EXPECT_NEAR(Kernel::strainSummedLogLikelihood(sums, 2, table), 0.0, 1e-12);
    EXPECT_NEAR(Kernel::strainSummedLogLikelihood(sums, 1, table), 0.0, 1e-12);
}

TEST(TransmissionKernelTest, LatentLocusTableMatchesThePlainComputation) {
    Engine rng(20261014);
    Kernel kernel;
    int compared = 0;
    for (int trial = 0; trial < 40; ++trial) {
        const unsigned int A = 2 + trial % 6;
        const auto freqs     = randomFrequencies(rng, A);
        const Genetics child = fromMask(A, randomMask(rng, A, 0.4));
        std::vector<Genetics> observed;
        for (unsigned int p = 0; p < static_cast<unsigned int>(trial % 3); ++p) {
            observed.push_back(fromMask(A, randomMask(rng, A, 0.4)));
        }
        const auto parents    = pointers(observed);
        const Variant variant = observed.empty() ? Variant::LatentParentOnly : Variant::ObservedAndLatentParents;
        const marginal::SourceOverlap overlap(freqs, Kernel::allelesOf(child), MAX_COI);

        const auto expected = marginal::locusTable<MAX_COI, MAX_STRAINS, Genetics>(overlap, A, [&](const Genetics& latent) -> std::optional<std::array<double, MAX_STRAINS>> {
            const auto term = kernel.locusTerm(variant, child, parents, &latent);
            if (term.impossible) {
                return std::nullopt;
            }
            return term.terms;
        });
        Kernel::LocusTable table;
        Kernel::contract(kernel.fTable(child, parents), overlap, table);

        const bool anyPossible = std::ranges::any_of(expected, [](const double x) { return x > NEG_INF; });
        EXPECT_EQ(table.impossible, not anyPossible) << "trial " << trial;
        if (table.impossible) {
            continue;
        }
        for (int k = 1; k <= MAX_COI; ++k) {
            for (int m = 1; m <= MAX_STRAINS; ++m) {
                const int cell = k * MAX_STRAINS + m - 1;
                if (k >= table.kMin and m >= table.mMin) {
                    EXPECT_NEAR(table.terms[cell], expected[cell], 1e-9 * std::max(1.0, std::abs(expected[cell]))) << "trial " << trial << " k=" << k << " m=" << m;
                    ++compared;
                } else {
                    EXPECT_EQ(expected[cell], NEG_INF) << "trial " << trial << " k=" << k << " m=" << m;
                }
            }
        }
    }
    EXPECT_GT(compared, 1000);
}

TEST(TransmissionKernelTest, LatentSumMatchesThePlainCombination) {
    Engine rng(7);
    Kernel kernel;
    const auto table = Kernel::strainCountTable(1.7);
    std::array<double, MAX_COI + 1> coiLogProbs{};
    coiLogProbs[0] = NEG_INF;
    double norm    = 0;
    for (int k = 1; k <= MAX_COI; ++k) {
        coiLogProbs[k] = -0.7 * k;
        norm += std::exp(coiLogProbs[k]);
    }
    for (int k = 1; k <= MAX_COI; ++k) {
        coiLogProbs[k] -= std::log(norm);
    }
    for (int trial = 0; trial < 20; ++trial) {
        const unsigned int A    = 3 + trial % 4;
        const Genetics child    = fromMask(A, randomMask(rng, A, 0.5));
        const Genetics observed = fromMask(A, randomMask(rng, A, 0.5));
        const std::vector<const Genetics*> parents{&observed};
        const marginal::SourceOverlap overlap(randomFrequencies(rng, A), Kernel::allelesOf(child), MAX_COI);
        Kernel::LocusTable locus;
        Kernel::contract(kernel.fTable(child, parents), overlap, locus);
        if (locus.impossible) {
            continue;
        }
        Kernel::LatentSums sums{};
        sums.add(locus);
        marginal::LocusTable<MAX_COI, MAX_STRAINS> plain{};
        plain.fill(NEG_INF);
        for (int k = locus.kMin; k <= MAX_COI; ++k) {
            for (int m = locus.mMin; m <= MAX_STRAINS; ++m) {
                plain[k * MAX_STRAINS + m - 1] = locus.terms[k * MAX_STRAINS + m - 1];
            }
        }
        const double expected = marginal::combine<MAX_COI, MAX_STRAINS>(plain, coiLogProbs, [&](const int m, const int n) { return Kernel::probNumStrains(table, m, n); }, 2, 0.0);
        EXPECT_NEAR(Kernel::latentSummedLogLikelihood(sums, 2, coiLogProbs, table), expected, 1e-9 * std::max(1.0, std::abs(expected))) << "trial " << trial;
    }
}
