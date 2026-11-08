//
// The family moves must leave their target invariant. Exact draws from a target over two infections' durations and
// genotypes are each moved by one GuidedGroupSampler relocation, or one PairLocusGibbsSampler step, and must still
// follow it. The harness is GuidedPairSamplerTest's.
//

#include "gtest/gtest.h"

#include "core/samplers/genetics/PairLocusGibbsSampler.h"
#include "core/samplers/specialized/Family.h"
#include "core/samplers/specialized/GuidedGroupSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"
#include "../../../impl/Studies.h"

#include <boost/random.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;
using transmission_nets::core::samplers::genetics::PairLocusGibbsSampler;
using transmission_nets::core::samplers::specialized::GuidedGroupSampler;
using transmission_nets::core::samplers::specialized::closePairs;
using transmission_nets::core::samplers::specialized::families;

namespace {
    using Engine = boost::random::mt19937;

    constexpr int DAYS                  = 60;
    constexpr unsigned int STRAIN_LIMIT = 2;

    // x (observed on day 400) and y (day 420) may infect each other; between them they were observed with all three
    // alleles of the locus.
    Study makeStudy() {
        Study study{.loci = {uniformLocus("L0", 3)}, .symptomaticInfectionDurations = flatDurations(DAYS), .asymptomaticInfectionDurations = flatDurations(DAYS)};
        study.infections.push_back({.id = "x", .observationTime = 400, .observedGenotypes = {{"L0", "110"}}, .permittedParents = {"y"}});
        study.infections.push_back({.id = "y", .observationTime = 420, .observedGenotypes = {{"L0", "011"}}, .permittedParents = {"x"}});
        return study;
    }

    // A duration's 10-day bin: bin j holds days 10j + 1 to 10j + 10.
    int binOf(const double duration) {
        return (static_cast<int>(std::round(duration)) - 1) / 10;
    }

    std::pair<double, double> binRange(const int j) {
        return {j == 0 ? 1.0 : 10.0 * j + 0.5, j == 5 ? DAYS - 0.5 : 10.0 * j + 10.5};
    }

    // x came first by bins when y's bin is less than two after x's (x was observed 20 days, two bins, earlier).
    bool xFirst(const int jx, const int jy) {
        return jy - jx < 2;
    }

    // The genotypes the target allows: one or two of the three alleles, as bit strings, and their index.
    const std::array<std::string, 6> GENOTYPES{"100", "010", "001", "110", "101", "011"};

    int genotypeIndex(const Genotype& g) {
        for (size_t k = 0; k < GENOTYPES.size(); ++k) {
            if (g == Genotype(GENOTYPES[k])) {
                return static_cast<int>(k);
            }
        }
        return -1;
    }

    // Whether every allele of `child` is in `parent`.
    bool within(const int child, const int parent) {
        const auto& c = GENOTYPES[child];
        const auto& p = GENOTYPES[parent];
        for (size_t a = 0; a < 3; ++a) {
            if (c[a] == '1' and p[a] == '0') {
                return false;
            }
        }
        return true;
    }

    // The part of the target coupling order and genotypes: whoever came first should carry the other's alleles.
    double coupling(const int jx, const int jy, const int gx, const int gy) {
        return xFirst(jx, jy) ? (within(gy, gx) ? 2.0 : -1.5) : (within(gx, gy) ? 1.5 : -2.0);
    }

    double logTarget(const double dx, const double dy, const int gx, const int gy) {
        if (dx < 1 or dy < 1 or dx >= DAYS - 0.5 or dy >= DAYS - 0.5 or gx < 0 or gy < 0) {
            return -std::numeric_limits<double>::infinity();
        }
        constexpr std::array<double, 6> wx{0.4, -0.3, 0.8, 0.1, -0.6, 0.2}, wy{-0.2, 0.5, 0.0, 0.7, -0.4, 0.3};
        constexpr std::array<double, 6> wg{0.3, -0.2, 0.1, 0.6, -0.4, 0.0};
        const int jx = binOf(dx), jy = binOf(dy);
        return wx[jx] + wy[jy] + wg[gx] - wg[gy] + coupling(jx, jy, gx, gy);
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261018);
        std::shared_ptr<State> state;
        std::shared_ptr<LocusImpl> locus;
        std::shared_ptr<InfectionEvent> x, y;

        Fixture() {
            state = std::make_shared<State>(makeStudy(), rng);
            locus = state->loci.at("L0");
            x     = state->infection("x");
            y     = state->infection("y");
        }

        [[nodiscard]] int gx() const { return genotypeIndex(x->latentGenotype(locus)->value()); }
        [[nodiscard]] int gy() const { return genotypeIndex(y->latentGenotype(locus)->value()); }
    };

    struct Target {
        const Fixture* f;
        double value() const { return logTarget(f->x->infectionDuration()->value(), f->y->infectionDuration()->value(), f->gx(), f->gy()); }
        [[nodiscard]] bool isDirty() const { return false; }
    };

    using Groups    = GuidedGroupSampler<Target, Engine, InfectionEvent, LocusImpl>;
    using PairLocus = PairLocusGibbsSampler<Target, Engine, InfectionEvent, LocusImpl>;

    struct Visits {
        double binDistance;
        double genotypeDistance;
        int reordered;
        int changedGenotypes;
    };

    // Exact draws from the target, each moved by one `update` and recorded by bins and by order and genotypes.
    Visits moveExactDraws(Fixture& f, const std::function<void()>& update, const int draws) {
        using Cell = std::tuple<int, int, int, int>;
        std::vector<Cell> cells;
        std::vector<double> cumulative;
        double total = 0;
        for (int jx = 0; jx < 6; ++jx) {
            for (int jy = 0; jy < 6; ++jy) {
                for (int gx = 0; gx < 6; ++gx) {
                    for (int gy = 0; gy < 6; ++gy) {
                        const auto [x0, x1] = binRange(jx);
                        const auto [y0, y1] = binRange(jy);
                        total += std::exp(logTarget((x0 + x1) / 2, (y0 + y1) / 2, gx, gy)) * (x1 - x0) * (y1 - y0);
                        cells.emplace_back(jx, jy, gx, gy);
                        cumulative.push_back(total);
                    }
                }
            }
        }
        using BinKey  = std::pair<int, int>;
        using GenoKey = std::tuple<bool, int, int>;
        std::map<BinKey, double> expectedBins;
        std::map<GenoKey, double> expectedGenotypes;
        for (size_t c = 0; c < cells.size(); ++c) {
            const auto [jx, jy, gx, gy] = cells[c];
            const double p              = (cumulative[c] - (c == 0 ? 0.0 : cumulative[c - 1])) / total;
            expectedBins[{jx, jy}] += p;
            expectedGenotypes[{xFirst(jx, jy), gx, gy}] += p;
        }

        boost::random::uniform_01<> uniform;
        std::map<BinKey, int> binVisits;
        std::map<GenoKey, int> genotypeVisits;
        Visits out{};
        for (int i = 0; i < draws; ++i) {
            const double u              = uniform(*f.rng) * total;
            const auto [jx, jy, gx, gy] = cells[std::min(cells.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
            const auto [x0, x1]         = binRange(jx);
            const auto [y0, y1]         = binRange(jy);
            f.x->infectionDuration()->initializeValue(x0 + uniform(*f.rng) * (x1 - x0));
            f.y->infectionDuration()->initializeValue(y0 + uniform(*f.rng) * (y1 - y0));
            f.x->latentGenotype(f.locus)->initializeValue(Genotype(GENOTYPES[gx]));
            f.y->latentGenotype(f.locus)->initializeValue(Genotype(GENOTYPES[gy]));

            update();

            const int kx = binOf(f.x->infectionDuration()->value()), ky = binOf(f.y->infectionDuration()->value());
            EXPECT_GE(f.gx(), 0);
            EXPECT_GE(f.gy(), 0);
            out.reordered += xFirst(kx, ky) != xFirst(jx, jy);
            out.changedGenotypes += f.gx() != gx or f.gy() != gy;
            ++binVisits[{kx, ky}];
            ++genotypeVisits[{xFirst(kx, ky), f.gx(), f.gy()}];
        }

        auto distance = [&](const auto& expected, const auto& visits) {
            double d = 0;
            for (const auto& [k, p] : expected) {
                d += std::abs((visits.contains(k) ? visits.at(k) : 0) / static_cast<double>(draws) - p) / 2;
            }
            return d;
        };
        out.binDistance      = distance(expectedBins, binVisits);
        out.genotypeDistance = distance(expectedGenotypes, genotypeVisits);
        return out;
    }

    std::vector<PairLocus::Pair> thePair(const Fixture& f) {
        return {{{f.x, DAYS - 0.5}, {f.y, DAYS - 0.5}}};
    }
}// namespace

TEST(FamilyMovesTest, GuidedGroupRelocationLeavesTargetInvariant) {
    Fixture f;
    auto target            = std::make_shared<Target>(Target{&f});
    const auto groups      = families(closePairs(thePair(f), 3));
    ASSERT_EQ(groups.size(), 1u);
    Groups sampler(groups, target, f.rng, 1.0, STRAIN_LIMIT, 10.0, 2, 1);
    constexpr int draws    = 300000;
    const auto visits      = moveExactDraws(f, [&] { sampler.update(); }, draws);
    EXPECT_GT(sampler.counts().accepted, sampler.counts().proposed / 20);
    EXPECT_GT(visits.reordered, draws / 50);
    EXPECT_LT(visits.binDistance, 0.01);
    EXPECT_LT(visits.genotypeDistance, 0.01);
}

TEST(FamilyMovesTest, PairLocusGibbsLeavesTargetInvariant) {
    Fixture f;
    auto target         = std::make_shared<Target>(Target{&f});
    PairLocus sampler(thePair(f), target, f.rng, STRAIN_LIMIT, 1, 3);
    ASSERT_EQ(sampler.units(), 1u);
    constexpr int draws = 300000;
    const auto visits   = moveExactDraws(f, [&] { sampler.update(); }, draws);
    EXPECT_GT(visits.changedGenotypes, draws / 20);
    EXPECT_LT(visits.binDistance, 0.01);
    EXPECT_LT(visits.genotypeDistance, 0.01);
}

// A Locus where either Infection has no Allele observed counts toward neither direction's novelty. x and y differ at
// L0 only; x is missing at L1 to L4 and y at L5 to L8, which once put both directions over the limit and left the
// pair without pair-locus units.
TEST(FamilyMovesTest, MissingLociDoNotCountTowardNovelty) {
    std::vector<std::shared_ptr<LocusImpl>> loci;
    auto x = std::make_shared<InfectionEvent>("x", 400, true);
    auto y = std::make_shared<InfectionEvent>("y", 420, true);
    for (int l = 0; l < 9; ++l) {
        loci.push_back(std::make_shared<LocusImpl>("L" + std::to_string(l), 3));
        const std::string ox = l == 0 ? "110" : (l <= 4 ? "000" : "100");
        const std::string oy = l == 0 ? "011" : (l <= 4 ? "100" : "000");
        x->addGenetics(loci.back(), Genotype(ox), Genotype(l == 0 ? ox : "100"));
        y->addGenetics(loci.back(), Genotype(oy), Genotype(l == 0 ? oy : "100"));
    }
    EXPECT_EQ(transmission_nets::core::samplers::specialized::novelty(*x, *y), 1u);
    EXPECT_EQ(transmission_nets::core::samplers::specialized::novelty(*y, *x), 1u);
    const std::vector<PairLocus::Pair> pairs{{{x, DAYS - 0.5}, {y, DAYS - 0.5}}};
    EXPECT_EQ(closePairs(pairs, 0).size(), 0u);
    EXPECT_EQ(closePairs(pairs, 1).size(), 1u);
    PairLocus sampler(pairs, nullptr, nullptr, STRAIN_LIMIT, 1, 3);
    EXPECT_EQ(sampler.units(), 1u);
}
