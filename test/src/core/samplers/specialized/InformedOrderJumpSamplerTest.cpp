//
// InformedOrderJumpSampler must leave its target invariant. y has three partners whose infection times split y's
// durations into four segments, each with its own weight, so the move has several tries to choose between. Exact draws
// from the target are each moved once and must still follow the target.
//

#include "gtest/gtest.h"

#include "core/samplers/specialized/InformedOrderJumpSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::samplers::specialized::InformedOrderJumpSampler;

namespace {
    using Engine = boost::random::mt19937;

    constexpr int DAYS = 60;
    // y is observed on day 450; partners a, b and c were infected on days 400, 415 and 430, so y came before a above a
    // duration of 50, before b above 35 and before c above 20.
    const std::vector<double> BOUNDARIES{20, 35, 50};

    double logTarget(const double yDuration) {
        double out = std::log(std::exp(-std::round(yDuration) / 25.0));// the duration's own weight, unlike the proposal's
        out += yDuration > 20 ? 1.2 : 0.0;
        out += yDuration > 35 ? -0.6 : 0.0;
        out += yDuration > 50 ? 0.9 : 0.0;
        return out;
    }

    int segment(const double d) {
        return static_cast<int>(std::ranges::count_if(BOUNDARIES, [&](const double b) { return d > b; }));
    }

    nlohmann::json makeInput() {
        nlohmann::json input;
        input["loci"].push_back({{"locus", "L0"}, {"num_alleles", 3}, {"allele_freqs", {1.0 / 3, 1.0 / 3, 1.0 / 3}}});
        auto node = [&](const std::string& id, const double time, const std::vector<std::string>& parents) {
            input["nodes"].push_back({{"id", id},
                                      {"observation_time", time},
                                      {"symptomatic", true},
                                      {"observed_genotype", {{{"locus", "L0"}, {"genotype", "110"}}}},
                                      {"allowed_parents", parents}});
        };
        node("a", 410, {"y"});
        node("b", 425, {"y"});
        node("c", 440, {"y"});
        node("y", 450, {"a", "b", "c"});
        return input;
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261012);
        std::shared_ptr<State> state;
        std::shared_ptr<InfectionEvent> y;
        std::vector<std::shared_ptr<InfectionEvent>> partners;

        Fixture() {
            const std::vector<double> idp(DAYS, 1.0 / DAYS);
            state = std::make_shared<State>(makeInput(), idp, idp, rng);
            for (const auto& id : {"a", "b", "c"}) {
                auto infection = state->infections[index(id)];
                infection->infectionDuration()->initializeValue(10);
                partners.push_back(infection);
            }
            y = state->infections[index("y")];
        }

        size_t index(const std::string& id) const {
            for (size_t i = 0; i < state->infections.size(); ++i) {
                if (state->infections[i]->id() == id) {
                    return i;
                }
            }
            throw std::out_of_range(id);
        }
    };

    struct Target {
        const Fixture* f;
        double value() const { return logTarget(f->y->infectionDuration()->value()); }
        [[nodiscard]] bool isDirty() const { return false; }
    };
}// namespace

// 0 tries every partner; 2 tries a random two of the three each update.
class InformedOrderJumpInvarianceTest : public testing::TestWithParam<size_t> {};

INSTANTIATE_TEST_SUITE_P(MaxTries, InformedOrderJumpInvarianceTest, testing::Values(size_t{0}, size_t{2}), [](const auto& info) { return std::to_string(info.param); });

TEST_P(InformedOrderJumpInvarianceTest, LeavesTargetInvariant) {
    Fixture f;
    auto target = std::make_shared<Target>(Target{&f});
    // Durations are proposed from a shape unlike the target, which the correction must undo.
    std::vector<double> shape(DAYS);
    for (size_t day = 0; day < shape.size(); ++day) {
        shape[day] = 1.0 + static_cast<double>(day);
    }
    InformedOrderJumpSampler<Target, Engine, InfectionEvent> sampler(f.y, f.partners, target, f.rng, 1.0, DAYS, shape, GetParam());

    // Exact draws: pieces of [1, 59.5) by weight (whole days, split at the boundaries), then uniformly within a piece.
    struct Cell {
        double first, second, weight;
    };
    std::vector<Cell> cells;
    std::vector<double> cumulative;
    double total = 0;
    for (int day = 1; day < DAYS; ++day) {
        std::vector<double> cuts{std::max(1.0, day - 0.5)};
        for (const double b : BOUNDARIES) {
            if (b > day - 0.5 and b < day + 0.5) {
                cuts.push_back(b);
            }
        }
        cuts.push_back(day + 0.5);
        for (size_t c = 0; c + 1 < cuts.size(); ++c) {
            const double first = cuts[c], second = cuts[c + 1];
            const double weight = std::exp(logTarget((first + second) / 2)) * (second - first);
            cells.push_back({first, second, weight});
            total += weight;
            cumulative.push_back(total);
        }
    }

    // States compared by segment and by day to within 5 days, so each piece lies within one bin.
    using Key = std::pair<int, int>;
    auto key = [&](const double d) { return Key{segment(d), static_cast<int>(std::round(d)) / 5}; };
    std::map<Key, double> expected;
    for (const auto& cell : cells) {
        expected[key((cell.first + cell.second) / 2)] += cell.weight / total;
    }

    boost::random::uniform_01<> uniform;
    constexpr int draws = 400000;
    std::map<Key, int> visits;
    int moved = 0, crossedSeveral = 0;
    for (int i = 0; i < draws; ++i) {
        const double u   = uniform(*f.rng) * total;
        const auto& cell = cells[std::min(cells.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        const double before = cell.first + uniform(*f.rng) * (cell.second - cell.first);
        f.y->infectionDuration()->initializeValue(before);
        sampler.update();
        const double after = f.y->infectionDuration()->value();
        if (after != before) {
            ++moved;
            ASSERT_NE(segment(after), segment(before));// every try crosses a partner
            crossedSeveral += std::abs(segment(after) - segment(before)) > 1;
        }
        ++visits[key(after)];
    }

    double distance = 0;
    for (const auto& [k, p] : expected) {
        distance += std::abs((visits.contains(k) ? visits.at(k) : 0) / static_cast<double>(draws) - p) / 2;
    }
    EXPECT_GT(moved, draws / 10);
    EXPECT_GT(crossedSeveral, draws / 100);// some moves jump past more than one partner
    EXPECT_LT(distance, 0.01);
}
