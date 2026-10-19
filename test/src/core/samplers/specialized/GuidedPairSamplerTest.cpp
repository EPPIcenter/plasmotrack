//
// GuidedPairSampler must leave its target invariant. Exact draws from a target over both infections' durations and
// genotypes are each moved by one guided swap and must still follow it.
//

#include "gtest/gtest.h"

#include "core/samplers/specialized/GuidedPairSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"
#include "../../../impl/Studies.h"

#include <boost/random.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;
using transmission_nets::core::samplers::specialized::GuidedPairSampler;

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

    using Sampler = GuidedPairSampler<Target, Engine, InfectionEvent, LocusImpl>;
}// namespace

TEST(GuidedPairSamplerTest, LeavesTargetInvariant) {
    Fixture f;
    auto target = std::make_shared<Target>(Target{&f});
    Sampler sampler(Sampler::mutualPairs(f.state->infections, *f.state->allowedRelationships, [](size_t) { return DAYS - 0.5; }), target, f.rng, 1.0, STRAIN_LIMIT);

    // Exact draws: two bins and two genotypes by weight, then each duration uniformly within its bin.
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
    // States are compared by bins, and by order and genotypes.
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
    constexpr int draws = 500000;
    std::map<BinKey, int> binVisits;
    std::map<GenoKey, int> genotypeVisits;
    int reordered = 0;
    for (int i = 0; i < draws; ++i) {
        const double u              = uniform(*f.rng) * total;
        const auto [jx, jy, gx, gy] = cells[std::min(cells.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        const auto [x0, x1]         = binRange(jx);
        const auto [y0, y1]         = binRange(jy);
        f.x->infectionDuration()->initializeValue(x0 + uniform(*f.rng) * (x1 - x0));
        f.y->infectionDuration()->initializeValue(y0 + uniform(*f.rng) * (y1 - y0));
        f.x->latentGenotype(f.locus)->initializeValue(Genotype(GENOTYPES[gx]));
        f.y->latentGenotype(f.locus)->initializeValue(Genotype(GENOTYPES[gy]));

        sampler.update();

        const int kx = binOf(f.x->infectionDuration()->value()), ky = binOf(f.y->infectionDuration()->value());
        ASSERT_GE(f.gx(), 0);
        ASSERT_GE(f.gy(), 0);
        reordered += xFirst(kx, ky) != xFirst(jx, jy);
        ++binVisits[{kx, ky}];
        ++genotypeVisits[{xFirst(kx, ky), f.gx(), f.gy()}];
    }

    auto distance = [&](const auto& expected, const auto& visits) {
        double out = 0;
        for (const auto& [k, p] : expected) {
            out += std::abs((visits.contains(k) ? visits.at(k) : 0) / static_cast<double>(draws) - p) / 2;
        }
        return out;
    };
    EXPECT_GT(sampler.counts().accepted, sampler.counts().proposed / 10);
    EXPECT_GT(reordered, draws / 20);
    EXPECT_LT(distance(expectedBins, binVisits), 0.01);
    EXPECT_LT(distance(expectedGenotypes, genotypeVisits), 0.01);
}

// Adapting sets the swaps an update makes from its evaluation budget, and a restore continues with the same pick
// probabilities, swaps per update and stats.
TEST(GuidedPairSamplerTest, RestoresAdaptation) {
    Fixture f;
    auto target         = std::make_shared<Target>(Target{&f});
    const auto pairs    = Sampler::mutualPairs(f.state->infections, *f.state->allowedRelationships, [](size_t) { return DAYS - 0.5; });
    const auto settings = Sampler::Settings{.moves = 1, .uniformShare = 0.5, .evaluations = 200};
    // x infected on day 370 and y on day 380, so their infection times can be exchanged.
    f.x->infectionDuration()->initializeValue(30);
    f.y->infectionDuration()->initializeValue(40);
    Sampler sampler(pairs, target, f.rng, 1.0, STRAIN_LIMIT, settings);
    EXPECT_EQ(sampler.adaptationState(), transmission_nets::core::samplers::AdaptationState{});
    for (int i = 0; i < 200; ++i) {
        sampler.update();
        sampler.adapt();
    }
    const auto [pair, different, stats, pick] = sampler.pairs().front();
    ASSERT_EQ(stats.swaps, sampler.counts().proposed);
    EXPECT_DOUBLE_EQ(pick, 1.0);
    EXPECT_GT(sampler.moves(), 1);
    const double n = static_cast<double>(stats.swaps);
    EXPECT_EQ(sampler.moves(), static_cast<int>(std::lround(200 / std::max(1.0, (stats.evaluations + 2 * stats.evaluations / n) / (n + 2)))));

    const auto state = sampler.adaptationState();
    EXPECT_EQ(state.updates, 200U);
    Sampler restored(pairs, target, f.rng, 1.0, STRAIN_LIMIT, settings);
    restored.restoreAdaptationState(state);
    EXPECT_EQ(restored.adaptationState(), state);
    EXPECT_EQ(restored.moves(), sampler.moves());
    const auto [rpair, rdifferent, rstats, rpick] = restored.pairs().front();
    EXPECT_EQ(rstats.swaps, stats.swaps);
    EXPECT_EQ(rstats.accepted, stats.accepted);
    EXPECT_DOUBLE_EQ(rpick, pick);

    auto wrong = state;
    wrong.fit.pop_back();
    EXPECT_THROW(restored.restoreAdaptationState(wrong), std::invalid_argument);
}
