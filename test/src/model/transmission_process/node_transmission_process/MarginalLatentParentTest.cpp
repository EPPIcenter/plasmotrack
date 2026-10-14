//
// Checks the exact expectation over a Source-drawn latent parent (MarginalLatentParentReference.h) against brute-force
// enumeration of every latent-parent genome on small loci.
//

#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/parameters/Parameter.h"

#include "MarginalLatentParentReference.h"
#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h"

#include <boost/random.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <vector>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

namespace {
    constexpr int MAX_ALLELES = 128;
    constexpr int MAX_COI     = 20;
    constexpr int MAX_PARENTS = 2;
    constexpr int MAX_STRAINS = 12;

    using Genetics   = AllelesBitSet<MAX_ALLELES>;
    using COIProb    = ZTPoisson<MAX_COI>;
    using PSP        = ZTGeometric<MAX_PARENTS + 1>;
    using Ntp        = MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, PSP, Genetics>;
    using Variant    = Ntp::Variant;
    using Engine     = boost::random::mt19937;
    using Table      = marginal::LocusTable<MAX_COI, MAX_STRAINS>;
    constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

    double logAdd(const double a, const double b) {
        const double m = std::max(a, b);
        return m == NEG_INF ? m : m + std::log(std::exp(a - m) + std::exp(b - m));
    }

    // P(distinct alleles of k draws = G) by inclusion-exclusion over subsets of G.
    double sourceGenotypeProb(const std::vector<double>& freqs, const unsigned int gMask, const int k) {
        double total = 0;
        for (unsigned int v = gMask;; v = (v - 1) & gMask) {
            double pv = 0;
            for (unsigned int i = 0; i < freqs.size(); ++i) {
                if (v >> i & 1u) {
                    pv += freqs[i];
                }
            }
            const int sign = (std::popcount(gMask) - std::popcount(v)) % 2 == 0 ? 1 : -1;
            total += sign * std::pow(pv, k);
            if (v == 0) {
                break;
            }
        }
        return total;
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

    std::vector<unsigned int> allelesOf(const unsigned int numAlleles, const unsigned int mask) {
        std::vector<unsigned int> out;
        for (unsigned int i = 0; i < numAlleles; ++i) {
            if (mask >> i & 1u) {
                out.push_back(i);
            }
        }
        return out;
    }
}// namespace

TEST(MarginalLatentParentTest, SourceOverlapMatchesEnumeration) {
    Engine rng(11);
    for (int trial = 0; trial < 20; ++trial) {
        const unsigned int A = 2 + trial % 6;
        const auto freqs     = randomFrequencies(rng, A);
        const unsigned int c = randomMask(rng, A, 0.5);
        const auto child     = allelesOf(A, c);
        const marginal::SourceOverlap overlap(freqs, child, MAX_COI);
        for (int k = 1; k <= MAX_COI; ++k) {
            // Aggregate P(G | k) by (U, w), with U as a mask over the child's alleles.
            std::vector<double> expected((1u << child.size()) * (A + 1), 0.0);
            double total = 0;
            for (unsigned int g = 1; g < (1u << A); ++g) {
                const double p = sourceGenotypeProb(freqs, g, k);
                total += p;
                unsigned int u = 0;
                for (unsigned int ci = 0; ci < child.size(); ++ci) {
                    if (g >> child[ci] & 1u) {
                        u |= 1u << ci;
                    }
                }
                const unsigned int w = std::popcount(g & ~c);
                expected[u * (A + 1) + w] += p;
            }
            EXPECT_NEAR(total, 1.0, 1e-9);
            for (unsigned int u = 1; u < (1u << child.size()); ++u) {
                for (unsigned int w = 0; w <= A - child.size(); ++w) {
                    const double got = w <= overlap.maxW() ? overlap(u, w, k) : 0.0;
                    EXPECT_NEAR(got, expected[u * (A + 1) + w], 1e-10) << "A=" << A << " k=" << k << " U=" << u << " w=" << w;
                }
            }
        }
    }
}

TEST(MarginalLatentParentTest, MatchesBruteForceOverLatentParentGenomes) {
    Engine rng(20261012);
    auto meanTransmitted = std::make_shared<Parameter<double>>(1.7);
    Ntp ntp(meanTransmitted);
    auto meanCoi = std::make_shared<Parameter<double>>(2.2);
    COIProb coiProb(meanCoi);
    auto psProb = std::make_shared<Parameter<double>>(0.6);
    PSP psp(psProb);
    const auto coiLogProbs = coiProb.value();

    double worst = 0;
    int cases    = 0;
    for (int trial = 0; trial < 60; ++trial) {
        const std::vector<unsigned int> alleleCounts = trial % 2 == 0 ? std::vector<unsigned int>{3, 4, 5} : std::vector<unsigned int>{2, 5, 4};
        const size_t L = alleleCounts.size();
        std::vector<std::vector<double>> freqs;
        std::vector<Genetics> child;
        std::vector<std::vector<Genetics>> parents(2);// [parent][locus]
        for (size_t l = 0; l < L; ++l) {
            const unsigned int A = alleleCounts[l];
            freqs.push_back(randomFrequencies(rng, A));
            child.push_back(fromMask(A, randomMask(rng, A, trial % 3 == 0 ? 0.6 : 0.3)));
            for (auto& p : parents) {
                p.push_back(fromMask(A, randomMask(rng, A, 0.4)));
            }
        }

        for (unsigned int numObserved = 0; numObserved <= 2; ++numObserved) {
            const int n           = static_cast<int>(numObserved) + 1;
            const Variant variant = numObserved == 0 ? Variant::LatentParentOnly : Variant::ObservedAndLatentParents;
            auto parentPtrs       = [&](const size_t l) {
                std::vector<const Genetics*> out;
                for (unsigned int p = 0; p < numObserved; ++p) {
                    out.push_back(&parents[p][l]);
                }
                return out;
            };

            // Exact.
            Table sums{};
            sums.fill(0.0);
            for (size_t l = 0; l < L; ++l) {
                const auto ptrs = parentPtrs(l);
                std::vector<unsigned int> childAlleles;
                for (unsigned int i = 0; i < alleleCounts[l]; ++i) {
                    if (child[l].allele(i)) {
                        childAlleles.push_back(i);
                    }
                }
                const marginal::SourceOverlap overlap(freqs[l], childAlleles, MAX_COI);
                const auto table = marginal::locusTable<MAX_COI, MAX_STRAINS, Genetics>(overlap, alleleCounts[l], [&](const Genetics& latent) -> std::optional<std::array<double, MAX_STRAINS>> {
                    const auto& term = ntp.locusTerm(variant, child[l], std::span<const Genetics* const>(ptrs), &latent);
                    if (term.impossible) {
                        return std::nullopt;
                    }
                    return term.terms;
                });
                marginal::accumulate<MAX_COI, MAX_STRAINS>(sums, table);
            }
            const double exact = marginal::combine<MAX_COI, MAX_STRAINS>(sums, coiLogProbs, [&](const int m, const int np) { return ntp.probNumStrains(m, np); }, n, std::log(psp.value()(n)));

            // Brute force: every latent-parent genome, its Source probability times the parent set's term.
            double brute = NEG_INF;
            std::vector<unsigned int> g(L, 1);
            while (true) {
                double source = NEG_INF;
                for (int k = 1; k <= MAX_COI; ++k) {
                    double lp = coiLogProbs[k];
                    for (size_t l = 0; l < L and lp > NEG_INF; ++l) {
                        const double p = sourceGenotypeProb(freqs[l], g[l], k);
                        lp += p > 0 ? std::log(p) : NEG_INF;
                    }
                    source = logAdd(source, lp);
                }
                Ntp::StrainSums strainSums{};
                for (size_t l = 0; l < L; ++l) {
                    const auto ptrs     = parentPtrs(l);
                    const Genetics lat  = fromMask(alleleCounts[l], g[l]);
                    strainSums.add(ntp.locusTerm(variant, child[l], std::span<const Genetics* const>(ptrs), &lat));
                }
                const double term = ntp.strainSummedLogLikelihood(strainSums, n);
                if (term > NEG_INF) {
                    brute = logAdd(brute, source + term + std::log(psp.value()(n)));
                }
                size_t l = 0;
                while (l < L and ++g[l] == (1u << alleleCounts[l])) {
                    g[l] = 1;
                    ++l;
                }
                if (l == L) {
                    break;
                }
            }

            ++cases;
            if (brute == NEG_INF) {
                EXPECT_EQ(exact, NEG_INF);
                continue;
            }
            EXPECT_NEAR(exact, brute, 1e-9) << "trial " << trial << " observed parents " << numObserved;
            worst = std::max(worst, std::abs(exact - brute));
        }
    }
    std::cout << "brute-force cases " << cases << ", worst |exact - brute| (log units) " << worst << "\n";
}
