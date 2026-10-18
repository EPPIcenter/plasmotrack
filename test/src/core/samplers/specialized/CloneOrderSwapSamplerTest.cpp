//
// CloneOrderSwapSampler must leave its target invariant. The move alone only alternates between mirrored states, so
// rather than running it as a chain, exact draws from the target are each moved once and must still follow the target.
//

#include "gtest/gtest.h"

#include "core/samplers/specialized/CloneOrderSwapSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"
#include "../../../impl/Studies.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;
using transmission_nets::core::parameters::Parameter;
using transmission_nets::core::samplers::specialized::CloneOrderSwapSampler;

namespace {
    using Engine = boost::random::mt19937;

    constexpr double Y_FIRST    = 0.9;// log weight of y having come first
    constexpr double X_DURATION = 10; // x was infected on day 390
    constexpr double BOUNDARY   = 30; // y (observed on day 420) came first above this
    constexpr int DAYS          = 60;

    // The duration's target density per day, unlike the proposal's.
    double durationWeight(const double duration) {
        return std::exp(-std::round(duration) / 20.0);
    }

    double logTarget(const double yDuration) {
        return (yDuration > BOUNDARY ? Y_FIRST : 0.0) + std::log(durationWeight(yDuration));
    }

    // One three-allele locus. x and y may infect each other; z may be infected by y but not infect it.
    Study makeStudy() {
        Study study{.loci = {uniformLocus("L0", 3)}, .symptomaticInfectionDurations = flatDurations(DAYS), .asymptomaticInfectionDurations = flatDurations(DAYS)};
        auto infection = [&](const std::string& id, const double time, const std::vector<std::string>& parents) {
            study.infections.push_back({.id = id, .observationTime = time, .observedGenotypes = {{"L0", "110"}}, .permittedParents = parents});
        };
        infection("x", 400, {"y"});
        infection("y", 420, {"x", "z"});
        infection("z", 410, {});
        return study;
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261009);
        std::shared_ptr<State> state;
        std::shared_ptr<LocusImpl> locus;
        std::shared_ptr<InfectionEvent> x, y;

        Fixture() {
            state   = std::make_shared<State>(makeStudy(), rng);
            locus   = state->loci.at("L0");
            x       = state->infection("x");
            y       = state->infection("y");
            x->infectionDuration()->initializeValue(X_DURATION);
            x->latentGenotype(locus)->initializeValue(Genotype("110"));
            y->latentGenotype(locus)->initializeValue(Genotype("110"));
        }
    };

    struct Target {
        const Fixture* f;

        double value() const {
            return logTarget(f->y->infectionDuration()->value());
        }

        [[nodiscard]] bool isDirty() const { return false; }
    };

    using Sampler = CloneOrderSwapSampler<Target, Engine, InfectionEvent>;

    // Moves y against its partners.
    Sampler makeSampler(const Fixture& f, std::shared_ptr<Target> target) {
        // Durations are proposed from a shape unlike the target, which the correction must undo.
        std::vector<double> shape(DAYS);
        for (size_t day = 0; day < shape.size(); ++day) {
            shape[day] = 1.0 + static_cast<double>(day);
        }
        return Sampler(f.y, Sampler::mutualPartners(f.state->infections, *f.state->allowedRelationships, f.state->indexOf("y")), std::move(target), f.rng, 1.0, DAYS, shape);
    }
}// namespace

TEST(CloneOrderSwapSamplerTest, PartnersCanInfectEachOther) {
    Fixture f;
    const auto partners = Sampler::mutualPartners(f.state->infections, *f.state->allowedRelationships, f.state->indexOf("y"));
    ASSERT_EQ(partners.size(), 1u);// not z, which cannot infect y
    EXPECT_EQ(partners[0].infection, f.x);
}

TEST(CloneOrderSwapSamplerTest, LeavesTargetInvariant) {
    Fixture f;
    auto target  = std::make_shared<Target>(Target{&f});
    auto sampler = makeSampler(f, target);

    // Exact draws: a piece of the duration's range by weight, then a duration uniformly within the piece. Durations
    // run over [1, 59.5), the days the proposal gives weight, split at the boundary.
    struct Cell {
        double first, second, weight;
    };
    std::vector<Cell> cells;
    std::vector<double> cumulative;
    double total = 0;
    for (int day = 1; day < DAYS; ++day) {
        const double low = std::max(1.0, day - 0.5), high = day + 0.5;
        for (const auto& [first, second] : {std::pair{low, std::min(high, BOUNDARY)}, std::pair{std::max(low, BOUNDARY), high}}) {
            if (second > first) {
                const double weight = std::exp(logTarget((first + second) / 2)) * (second - first);
                cells.push_back({first, second, weight});
                total += weight;
                cumulative.push_back(total);
            }
        }
    }

    // States are compared by which came first, and by which came first and the duration's day to within 10 days (so
    // each piece lies within one bin).
    using SideKey     = bool;
    using DurationKey = std::pair<bool, int>;
    auto sideKey      = [&] { return SideKey{f.y->infectionDuration()->value() > BOUNDARY}; };
    auto durationKey  = [&] {
        const double duration = f.y->infectionDuration()->value();
        return DurationKey{duration > BOUNDARY, static_cast<int>(std::round(duration)) / 10};
    };
    std::map<SideKey, double> expectedSides;
    std::map<DurationKey, double> expectedDurations;
    for (const auto& cell : cells) {
        f.y->infectionDuration()->initializeValue((cell.first + cell.second) / 2);
        expectedSides[sideKey()] += cell.weight / total;
        expectedDurations[durationKey()] += cell.weight / total;
    }

    boost::random::uniform_01<> uniform;
    constexpr int draws = 3000000;
    std::map<SideKey, int> sideVisits;
    std::map<DurationKey, int> durationVisits;
    int moved = 0;
    for (int i = 0; i < draws; ++i) {
        const double u   = uniform(*f.rng) * total;
        const auto& cell = cells[std::min(cells.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        const double before = cell.first + uniform(*f.rng) * (cell.second - cell.first);
        f.y->infectionDuration()->initializeValue(before);

        sampler.update();

        const double after = f.y->infectionDuration()->value();
        if (after != before) {
            // A move puts y on the other side of x; genotypes are untouched.
            ++moved;
            ASSERT_NE(after > BOUNDARY, before > BOUNDARY);
        }
        ASSERT_EQ(f.x->latentGenotype(f.locus)->value(), Genotype("110"));
        ASSERT_EQ(f.y->latentGenotype(f.locus)->value(), Genotype("110"));
        ++sideVisits[sideKey()];
        ++durationVisits[durationKey()];
    }

    auto distance = [&](const auto& expected, const auto& visits) {
        double out = 0;
        for (const auto& [k, p] : expected) {
            out += std::abs((visits.contains(k) ? visits.at(k) : 0) / static_cast<double>(draws) - p) / 2;
        }
        return out;
    };
    EXPECT_GT(moved, draws / 10);
    EXPECT_LT(distance(expectedSides, sideVisits), 0.005);
    EXPECT_LT(distance(expectedDurations, durationVisits), 0.005);
}

TEST(CloneOrderSwapSamplerTest, OnlyMovesAgainstClones) {
    Fixture f;
    f.x->latentGenotype(f.locus)->initializeValue(Genotype("100"));
    auto target  = std::make_shared<Target>(Target{&f});
    auto sampler = makeSampler(f, target);
    for (int i = 0; i < 100; ++i) {
        sampler.update();
    }
    EXPECT_EQ(sampler.counts().proposed, 0u);
}

namespace {
    // A target over both durations, in 10-day bins (bin j holds days 10j + 1 to 10j + 10): a weight for each one's
    // bin and one for the difference of the bins. Exchanging infection times (x observed 20 days before y) shifts each
    // duration by 20 days, two bins.
    double logPairTarget(const double xDuration, const double yDuration) {
        if (xDuration < 1 or yDuration < 1 or xDuration >= DAYS - 0.5 or yDuration >= DAYS - 0.5) {
            return -std::numeric_limits<double>::infinity();
        }
        constexpr std::array<double, 6> wx{0.4, -0.3, 0.8, 0.1, -0.6, 0.2}, wy{-0.2, 0.5, 0.0, 0.7, -0.4, 0.3};
        const int jx = (static_cast<int>(std::round(xDuration)) - 1) / 10, jy = (static_cast<int>(std::round(yDuration)) - 1) / 10;
        return wx[jx] + wy[jy] - 0.3 * std::abs(jy - jx - 2);
    }

    struct PairTarget {
        const Fixture* f;
        double value() const { return logPairTarget(f->x->infectionDuration()->value(), f->y->infectionDuration()->value()); }
        [[nodiscard]] bool isDirty() const { return false; }
    };
}// namespace

// With time exchange, the move still leaves a target over both clones' durations invariant, and exchanges happen.
TEST(CloneOrderSwapSamplerTest, TimeExchangeLeavesPairTargetInvariant) {
    Fixture f;
    auto target = std::make_shared<PairTarget>(PairTarget{&f});
    std::vector<double> shape(DAYS);
    for (size_t day = 0; day < shape.size(); ++day) {
        shape[day] = 1.0 + static_cast<double>(day);
    }
    auto partners = CloneOrderSwapSampler<PairTarget, Engine, InfectionEvent>::mutualPartners(f.state->infections, *f.state->allowedRelationships, f.state->indexOf("y"));
    for (auto& partner : partners) {
        partner.maxDuration = DAYS;
    }
    CloneOrderSwapSampler<PairTarget, Engine, InfectionEvent> sampler(f.y, std::move(partners), target, f.rng, 1.0, DAYS, shape);
    sampler.enableTimeExchange();

    // Exact draws: a pair of bins by weight, then each duration uniformly within its bin.
    auto range = [](const int j) { return std::pair{j == 0 ? 1.0 : 10.0 * j + 0.5, j == 5 ? DAYS - 0.5 : 10.0 * j + 10.5}; };
    std::vector<std::pair<int, int>> bins;
    std::vector<double> cumulative, probability;
    double total = 0;
    for (int jx = 0; jx < 6; ++jx) {
        for (int jy = 0; jy < 6; ++jy) {
            const auto [x0, x1] = range(jx);
            const auto [y0, y1] = range(jy);
            const double weight = std::exp(logPairTarget((x0 + x1) / 2, (y0 + y1) / 2)) * (x1 - x0) * (y1 - y0);
            bins.emplace_back(jx, jy);
            probability.push_back(weight);
            total += weight;
            cumulative.push_back(total);
        }
    }
    auto binOf = [](const double d) { return (static_cast<int>(std::round(d)) - 1) / 10; };

    boost::random::uniform_01<> uniform;
    constexpr int draws = 1000000;
    std::map<std::pair<int, int>, int> visits;
    int exchanged = 0;
    for (int i = 0; i < draws; ++i) {
        const double u     = uniform(*f.rng) * total;
        const auto [jx, jy] = bins[std::min(bins.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        const auto [x0, x1] = range(jx);
        const auto [y0, y1] = range(jy);
        const double xBefore = x0 + uniform(*f.rng) * (x1 - x0);
        f.x->infectionDuration()->initializeValue(xBefore);
        f.y->infectionDuration()->initializeValue(y0 + uniform(*f.rng) * (y1 - y0));

        sampler.update();

        exchanged += f.x->infectionDuration()->value() != xBefore;
        ++visits[{binOf(f.x->infectionDuration()->value()), binOf(f.y->infectionDuration()->value())}];
    }
    double distance = 0;
    for (size_t b = 0; b < bins.size(); ++b) {
        distance += std::abs((visits.contains(bins[b]) ? visits.at(bins[b]) : 0) / static_cast<double>(draws) - probability[b] / total) / 2;
    }
    EXPECT_GT(exchanged, draws / 10);
    EXPECT_LT(distance, 0.005);
}
