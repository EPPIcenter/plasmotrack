//
// Checks the exact expectation over a Source-drawn latent parent (MarginalLatentParentReference.h) against brute-force
// enumeration of every latent-parent genome on small loci.
//

#include "gtest/gtest.h"

#include "core/datatypes/Genotype.h"
#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/parameters/Parameter.h"

#include "MarginalLatentParentReference.h"
#include "model/transmission_process/TransmissionKernel.h"
#include <boost/random.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

namespace {
    constexpr int MAX_COI     = 20;
    constexpr int MAX_PARENTS = 2;
    constexpr int MAX_STRAINS = 12;
    using COIProb    = ZTPoisson<MAX_COI>;
    using PSP        = ZTGeometric<MAX_PARENTS + 1>;
    using Kernel     = TransmissionKernel<MAX_PARENTS, MAX_STRAINS, MAX_COI>;
    using Variant    = Kernel::Variant;
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

    Genotype fromMask(const unsigned int numAlleles, const unsigned int mask) {
        Genotype g(static_cast<int>(numAlleles));
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
            // U = ∅ is read with a transmission leak, by the latent parent alone; it needs w >= 1.
            for (unsigned int u = 0; u < (1u << child.size()); ++u) {
                for (unsigned int w = u == 0 ? 1 : 0; w <= A - child.size(); ++w) {
                    const double got = w <= overlap.maxW() ? overlap(u, w, k) : 0.0;
                    EXPECT_NEAR(got, expected[u * (A + 1) + w], 1e-10) << "A=" << A << " k=" << k << " U=" << u << " w=" << w;
                }
            }
        }
    }
}

TEST(MarginalLatentParentTest, MatchesBruteForceOverLatentParentGenomes) {
    Engine rng(20261012);
    Kernel kernel;
    const auto strainTable = Kernel::strainCountTable(1.7);
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
        std::vector<Genotype> child;
        std::vector<std::vector<Genotype>> parents(2);// [parent][locus]
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
                std::vector<const Genotype*> out;
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
                const auto table = marginal::locusTable<MAX_COI, MAX_STRAINS>(overlap, alleleCounts[l], [&](const Genotype& latent) -> std::optional<std::array<double, MAX_STRAINS>> {
                    const auto& term = kernel.locusTerm(variant, child[l], std::span<const Genotype* const>(ptrs), &latent);
                    if (term.impossible) {
                        return std::nullopt;
                    }
                    return term.terms;
                });
                marginal::accumulate<MAX_COI, MAX_STRAINS>(sums, table);
            }
            const double exact = marginal::combine<MAX_COI, MAX_STRAINS>(sums, coiLogProbs, [&](const int m, const int np) { return Kernel::probNumStrains(strainTable, m, np); }, n, std::log(psp.value()(n)));

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
                Kernel::StrainSums strainSums{};
                for (size_t l = 0; l < L; ++l) {
                    const auto ptrs     = parentPtrs(l);
                    const Genotype lat  = fromMask(alleleCounts[l], g[l]);
                    strainSums.add(kernel.locusTerm(variant, child[l], std::span<const Genotype* const>(ptrs), &lat));
                }
                const double term = Kernel::strainSummedLogLikelihood(strainSums, n, strainTable);
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

// The transmission leak ε > 0: the parent pool at a locus is (1 - ε) p + ε π. The kernel computes f by a generating
// function; the reference (leakyLocusLogTerm) follows the strains drawn one at a time.
namespace {
    constexpr double LEAK = 0.05;

    std::optional<std::array<double, MAX_STRAINS>> kernelTerm(Kernel& kernel, const Variant variant, const Genotype& child, const std::vector<const Genotype*>& observed, const Genotype* latent, const std::vector<double>& freqs) {
        const auto term = kernel.locusTerm(variant, child, std::span<const Genotype* const>(observed), latent, &freqs);
        if (term.impossible) {
            return std::nullopt;
        }
        return term.terms;
    }

    std::optional<std::array<double, MAX_STRAINS>> referenceTerm(const Variant variant, const Genotype& child, const std::vector<const Genotype*>& observed, const Genotype* latent, const std::vector<double>& freqs, const double leak) {
        auto parents = observed;
        if (latent != nullptr) {
            parents.push_back(latent);
        }
        return marginal::leakyLocusLogTerm<MAX_STRAINS>(child, parents, freqs, leak, variant == Variant::ObservedAndLatentParents);
    }

    void expectTermsNear(const std::optional<std::array<double, MAX_STRAINS>>& got, const std::optional<std::array<double, MAX_STRAINS>>& expected, const std::string& what) {
        ASSERT_EQ(got.has_value(), expected.has_value()) << what;
        if (not got) {
            return;
        }
        for (int m = 0; m < MAX_STRAINS; ++m) {
            if (std::isfinite((*expected)[m])) {
                EXPECT_NEAR((*got)[m], (*expected)[m], 1e-10 * std::max(1.0, std::abs((*expected)[m]))) << what << ", m = " << m + 1;
            } else {
                EXPECT_EQ((*got)[m], (*expected)[m]) << what << ", m = " << m + 1;
            }
        }
    }
}// namespace

TEST(MarginalLatentParentTest, LeakyLocusTermMatchesStrainByStrainDraws) {
    Engine rng(20261020);
    Kernel kernel(LEAK);
    for (int trial = 0; trial < 300; ++trial) {
        const unsigned int A = 2 + trial % 6;
        const auto freqs     = randomFrequencies(rng, A);
        const Genotype child = fromMask(A, randomMask(rng, A, 0.5));
        const Genotype p0    = fromMask(A, randomMask(rng, A, 0.4));
        const Genotype p1    = fromMask(A, randomMask(rng, A, 0.4));
        const Genotype lat   = fromMask(A, randomMask(rng, A, 0.5));
        const std::vector<const Genotype*> none{}, one{&p0}, two{&p0, &p1};
        const auto what = "trial " + std::to_string(trial);
        expectTermsNear(kernelTerm(kernel, Variant::ObservedParents, child, one, nullptr, freqs), referenceTerm(Variant::ObservedParents, child, one, nullptr, freqs, LEAK), what + ", one parent");
        expectTermsNear(kernelTerm(kernel, Variant::ObservedParents, child, two, nullptr, freqs), referenceTerm(Variant::ObservedParents, child, two, nullptr, freqs, LEAK), what + ", two parents");
        expectTermsNear(kernelTerm(kernel, Variant::ObservedAndLatentParents, child, one, &lat, freqs), referenceTerm(Variant::ObservedAndLatentParents, child, one, &lat, freqs, LEAK), what + ", one parent and L");
        expectTermsNear(kernelTerm(kernel, Variant::ObservedAndLatentParents, child, two, &lat, freqs), referenceTerm(Variant::ObservedAndLatentParents, child, two, &lat, freqs, LEAK), what + ", two parents and L");
        expectTermsNear(kernelTerm(kernel, Variant::LatentParentOnly, child, none, &lat, freqs), referenceTerm(Variant::LatentParentOnly, child, none, &lat, freqs, LEAK), what + ", L alone");
    }
}

// A child allele no parent carries is impossible without a leak and costs a finite amount with one. The rule that the
// latent parent shares an allele with the child, when with observed parents, holds either way.
TEST(MarginalLatentParentTest, ALeakMakesAnUncarriedChildAlleleFinite) {
    Kernel exact;
    Kernel leaky(LEAK);
    const std::vector<double> freqs{0.4, 0.3, 0.2, 0.1};
    const Genotype child  = fromMask(4, 0b0011);// alleles 0 and 1
    const Genotype parent = fromMask(4, 0b0001);// lacks allele 1
    const Genotype other  = fromMask(4, 0b1100);// shares nothing with the child
    const std::vector<const Genotype*> observed{&parent};
    const std::vector<const Genotype*> none{};

    EXPECT_TRUE(exact.locusTerm(Variant::ObservedParents, child, observed, nullptr, &freqs).impossible);
    const auto observedOnly = leaky.locusTerm(Variant::ObservedParents, child, observed, nullptr, &freqs);
    ASSERT_FALSE(observedOnly.impossible);
    EXPECT_TRUE(std::isfinite(observedOnly.terms[1]));
    expectTermsNear(observedOnly.terms, referenceTerm(Variant::ObservedParents, child, observed, nullptr, freqs, LEAK), "observed parent");

    // The latent parent alone, sharing nothing with the child.
    EXPECT_TRUE(exact.locusTerm(Variant::LatentParentOnly, child, none, &other, &freqs).impossible);
    const auto latentOnly = leaky.locusTerm(Variant::LatentParentOnly, child, none, &other, &freqs);
    ASSERT_FALSE(latentOnly.impossible);
    expectTermsNear(latentOnly.terms, referenceTerm(Variant::LatentParentOnly, child, none, &other, freqs, LEAK), "L alone");

    // With an observed parent, a latent parent sharing nothing with the child stays impossible.
    EXPECT_TRUE(leaky.locusTerm(Variant::ObservedAndLatentParents, child, observed, &other, &freqs).impossible);
    // One that shares allele 0 with the child, while neither carries allele 1, is finite only with the leak.
    const Genotype sharing = fromMask(4, 0b0101);
    EXPECT_TRUE(exact.locusTerm(Variant::ObservedAndLatentParents, child, observed, &sharing, &freqs).impossible);
    EXPECT_FALSE(leaky.locusTerm(Variant::ObservedAndLatentParents, child, observed, &sharing, &freqs).impossible);
}

// The model's path with a leak, F table and its contraction with the Source overlap summed over loci, against brute
// force over every latent-parent genome with the plain reference term. Also the observed parents alone, against the
// plain term. Parents and frequencies are drawn so that children often carry alleles no observed parent does.
TEST(MarginalLatentParentTest, LeakyMatchesBruteForceOverLatentParentGenomes) {
    Engine rng(20261021);
    Kernel kernel(LEAK);
    const auto strainTable = Kernel::strainCountTable(1.7);
    auto meanCoi           = std::make_shared<Parameter<double>>(2.2);
    COIProb coiProb(meanCoi);
    auto psProb = std::make_shared<Parameter<double>>(0.6);
    PSP psp(psProb);
    const auto coiLogProbs = coiProb.value();

    double worst      = 0;
    int cases         = 0;
    int uncarriedSets = 0;
    for (int trial = 0; trial < 40; ++trial) {
        const std::vector<unsigned int> alleleCounts = trial % 2 == 0 ? std::vector<unsigned int>{3, 4, 5} : std::vector<unsigned int>{2, 5, 4};
        const size_t L = alleleCounts.size();
        std::vector<std::vector<double>> freqs;
        std::vector<Genotype> child;
        std::vector<std::vector<Genotype>> parents(2);
        for (size_t l = 0; l < L; ++l) {
            const unsigned int A = alleleCounts[l];
            freqs.push_back(randomFrequencies(rng, A));
            child.push_back(fromMask(A, randomMask(rng, A, trial % 3 == 0 ? 0.6 : 0.3)));
            for (auto& p : parents) {
                p.push_back(fromMask(A, randomMask(rng, A, 0.3)));
            }
        }

        for (unsigned int numObserved = 0; numObserved <= 2; ++numObserved) {
            const int nWith       = static_cast<int>(numObserved) + 1;
            const Variant variant = numObserved == 0 ? Variant::LatentParentOnly : Variant::ObservedAndLatentParents;
            auto parentPtrs       = [&](const size_t l) {
                std::vector<const Genotype*> out;
                for (unsigned int p = 0; p < numObserved; ++p) {
                    out.push_back(&parents[p][l]);
                }
                return out;
            };

            // The model's path: F table, contraction, sums over loci, then the sum over k and m.
            Kernel::LatentSums latentSums{};
            for (size_t l = 0; l < L; ++l) {
                const auto ptrs         = parentPtrs(l);
                const auto childAlleles = Kernel::allelesOf(child[l]);
                const auto f            = kernel.fTable(child[l], std::span<const Genotype* const>(ptrs), &freqs[l]);
                Kernel::LocusTable table;
                Kernel::contract(f, marginal::SourceOverlap(freqs[l], childAlleles, MAX_COI), table);
                latentSums.add(table);
            }
            double exact = Kernel::latentSummedLogLikelihood(latentSums, nWith, coiLogProbs, strainTable);
            if (exact > NEG_INF) {
                exact += std::log(psp.value()(nWith));
            }

            // Brute force: every latent-parent genome, its Source probability times the plain term.
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
                Kernel::StrainSums strainSums{};
                for (size_t l = 0; l < L; ++l) {
                    const Genotype lat = fromMask(alleleCounts[l], g[l]);
                    const auto term    = referenceTerm(variant, child[l], parentPtrs(l), &lat, freqs[l], LEAK);
                    Kernel::LocusTerm locus{};
                    locus.impossible = not term.has_value();
                    if (term) {
                        locus.terms = *term;
                    }
                    strainSums.add(locus);
                }
                const double term = Kernel::strainSummedLogLikelihood(strainSums, nWith, strainTable);
                if (term > NEG_INF) {
                    brute = logAdd(brute, source + term + std::log(psp.value()(nWith)));
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
            } else {
                EXPECT_NEAR(exact, brute, 1e-9) << "trial " << trial << " observed parents " << numObserved << " with L";
                worst = std::max(worst, std::abs(exact - brute));
            }

            if (numObserved == 0) {
                continue;
            }
            // The observed parents alone: the kernel's terms summed over loci against the plain ones.
            Kernel::StrainSums modelSums{}, plainSums{};
            bool uncarried = false;
            for (size_t l = 0; l < L; ++l) {
                const auto ptrs = parentPtrs(l);
                modelSums.add(kernel.locusTerm(Variant::ObservedParents, child[l], std::span<const Genotype* const>(ptrs), nullptr, &freqs[l]));
                const auto term = referenceTerm(Variant::ObservedParents, child[l], ptrs, nullptr, freqs[l], LEAK);
                Kernel::LocusTerm locus{};
                locus.impossible = not term.has_value();
                if (term) {
                    locus.terms = *term;
                }
                plainSums.add(locus);
                for (unsigned int a = 0; a < alleleCounts[l]; ++a) {
                    uncarried = uncarried or (child[l].allele(a) and std::ranges::none_of(ptrs, [&](const Genotype* p) { return p->allele(a); }));
                }
            }
            const double model = Kernel::strainSummedLogLikelihood(modelSums, numObserved, strainTable);
            const double plain = Kernel::strainSummedLogLikelihood(plainSums, numObserved, strainTable);
            ++cases;
            if (uncarried) {
                ++uncarriedSets;
                EXPECT_TRUE(std::isfinite(model)) << "trial " << trial << ": a child allele no parent carries is finite with a leak";
            }
            if (std::isfinite(plain)) {
                EXPECT_NEAR(model, plain, 1e-9 * std::max(1.0, std::abs(plain))) << "trial " << trial << " observed parents " << numObserved;
                worst = std::max(worst, std::abs(model - plain));
            } else {
                EXPECT_EQ(model, plain);
            }
        }
    }
    EXPECT_GT(uncarriedSets, 10);
    std::cout << "leaky brute-force cases " << cases << " (" << uncarriedSets << " observed-only sets with an uncarried child allele), worst |exact - brute| (log units) " << worst << "\n";
}
