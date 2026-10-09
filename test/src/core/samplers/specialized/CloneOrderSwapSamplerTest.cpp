//
// CloneOrderSwapSampler must leave its target invariant. The move alone only alternates between mirrored states, so
// rather than running it as a chain, exact draws from the target are each moved once and must still follow the target.
//

#include "gtest/gtest.h"

#include "core/samplers/specialized/CloneOrderSwapSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::parameters::Parameter;
using transmission_nets::core::samplers::specialized::CloneOrderSwapSampler;

namespace {
    using Engine = boost::random::mt19937;

    constexpr std::array<double, 3> X_WEIGHTS{0.7, -0.4, 0.2};    // on x's latent parent
    constexpr std::array<double, 3> Y_WEIGHTS{-0.3, 0.5, 0.1};    // on y's latent parent
    constexpr std::array<double, 3> FIRST_WEIGHTS{1.0, 0.8, -0.6};// on the latent parent of whichever came first
    constexpr double X_DURATION = 10;                              // x was infected on day 390
    constexpr double BOUNDARY   = 30;                              // y (observed on day 420) came first above this
    constexpr int DAYS          = 60;

    const std::vector<std::string> GENOTYPES{"100", "010", "001", "110", "101", "011", "111"};

    double weigh(const std::array<double, 3>& weights, const GeneticsImpl& genotype) {
        double out = 0;
        for (unsigned int j = 0; j < weights.size(); ++j) {
            out += genotype.allele(j) ? weights[j] : 0.0;
        }
        return out;
    }

    // The duration's target density per day, unlike the proposal's.
    double durationWeight(const double duration) {
        return std::exp(-std::round(duration) / 20.0);
    }

    double logTarget(const GeneticsImpl& xLatent, const GeneticsImpl& yLatent, const double yDuration) {
        const bool yFirst = yDuration > BOUNDARY;
        return weigh(X_WEIGHTS, xLatent) + weigh(Y_WEIGHTS, yLatent) + weigh(FIRST_WEIGHTS, yFirst ? yLatent : xLatent) + std::log(durationWeight(yDuration));
    }

    // One three-allele locus. x and y may infect each other; z may be infected by y but not infect it.
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
        node("x", 400, {"y"});
        node("y", 420, {"x", "z"});
        node("z", 410, {});
        return input;
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261009);
        std::shared_ptr<State> state;
        std::shared_ptr<LocusImpl> locus;
        std::shared_ptr<InfectionEvent> x, y, xLatent, yLatent;

        Fixture() {
            const std::vector<double> idp(DAYS, 1.0 / DAYS);
            state   = std::make_shared<State>(makeInput(), idp, idp, rng);
            locus   = state->loci.at("L0");
            x       = state->infections[index("x")];
            y       = state->infections[index("y")];
            xLatent = state->latentParents[index("x")];
            yLatent = state->latentParents[index("y")];
            x->infectionDuration()->initializeValue(X_DURATION);
            x->latentGenotype(locus)->initializeValue(GeneticsImpl("110"));
            y->latentGenotype(locus)->initializeValue(GeneticsImpl("110"));
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

        double value() const {
            return logTarget(f->xLatent->latentGenotype(f->locus)->value(), f->yLatent->latentGenotype(f->locus)->value(), f->y->infectionDuration()->value());
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
        return Sampler(f.y, f.yLatent, Sampler::mutualPartners(f.state->infections, f.state->latentParents, *f.state->allowedRelationships, f.index("y")), std::move(target), f.rng, 1.0, DAYS, shape);
    }
}// namespace

TEST(CloneOrderSwapSamplerTest, PartnersCanInfectEachOther) {
    Fixture f;
    const auto partners = Sampler::mutualPartners(f.state->infections, f.state->latentParents, *f.state->allowedRelationships, f.index("y"));
    ASSERT_EQ(partners.size(), 1u);// not z, which cannot infect y
    EXPECT_EQ(partners[0].infection, f.x);
    EXPECT_EQ(partners[0].latentParent, f.xLatent);
}

TEST(CloneOrderSwapSamplerTest, LeavesTargetInvariant) {
    Fixture f;
    auto target  = std::make_shared<Target>(Target{&f});
    auto sampler = makeSampler(f, target);

    // Exact draws: latent parents' genotypes and a piece of the duration's range by weight, then a duration uniformly
    // within the piece. Durations run over [1, 59.5), the days the proposal gives weight, split at the boundary.
    struct Cell {
        GeneticsImpl xLatent, yLatent;
        double first, second, weight;
    };
    std::vector<Cell> cells;
    std::vector<double> cumulative;
    double total = 0;
    for (const auto& xg : GENOTYPES) {
        for (const auto& yg : GENOTYPES) {
            for (int day = 1; day < DAYS; ++day) {
                const double low = std::max(1.0, day - 0.5), high = day + 0.5;
                for (const auto& [first, second] : {std::pair{low, std::min(high, BOUNDARY)}, std::pair{std::max(low, BOUNDARY), high}}) {
                    if (second > first) {
                        const double mid    = (first + second) / 2;
                        const double weight = std::exp(logTarget(GeneticsImpl(xg), GeneticsImpl(yg), mid)) * (second - first);
                        cells.push_back({GeneticsImpl(xg), GeneticsImpl(yg), first, second, weight});
                        total += weight;
                        cumulative.push_back(total);
                    }
                }
            }
        }
    }

    // States are compared by their latent parents and which came first, and separately by which came first and the
    // duration's day to within 10 days (so each piece lies within one bin).
    using GeneticsKey = std::tuple<std::string, std::string, bool>;
    using DurationKey = std::pair<bool, int>;
    auto geneticsKey  = [&] {
        return GeneticsKey{f.xLatent->latentGenotype(f.locus)->value().serialize(), f.yLatent->latentGenotype(f.locus)->value().serialize(), f.y->infectionDuration()->value() > BOUNDARY};
    };
    auto durationKey = [&] {
        const double duration = f.y->infectionDuration()->value();
        return DurationKey{duration > BOUNDARY, static_cast<int>(std::round(duration)) / 10};
    };
    std::map<GeneticsKey, double> expectedGenetics;
    std::map<DurationKey, double> expectedDurations;
    for (const auto& cell : cells) {
        f.xLatent->latentGenotype(f.locus)->initializeValue(cell.xLatent);
        f.yLatent->latentGenotype(f.locus)->initializeValue(cell.yLatent);
        f.y->infectionDuration()->initializeValue((cell.first + cell.second) / 2);
        expectedGenetics[geneticsKey()] += cell.weight / total;
        expectedDurations[durationKey()] += cell.weight / total;
    }

    boost::random::uniform_01<> uniform;
    constexpr int draws = 3000000;
    std::map<GeneticsKey, int> geneticsVisits;
    std::map<DurationKey, int> durationVisits;
    int moved = 0;
    for (int i = 0; i < draws; ++i) {
        const double u   = uniform(*f.rng) * total;
        const auto& cell = cells[std::min(cells.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        f.xLatent->latentGenotype(f.locus)->initializeValue(cell.xLatent);
        f.yLatent->latentGenotype(f.locus)->initializeValue(cell.yLatent);
        const double before = cell.first + uniform(*f.rng) * (cell.second - cell.first);
        f.y->infectionDuration()->initializeValue(before);

        sampler.update();

        const double after = f.y->infectionDuration()->value();
        const auto& xNow   = f.xLatent->latentGenotype(f.locus)->value();
        const auto& yNow   = f.yLatent->latentGenotype(f.locus)->value();
        if (after != before) {
            // A move puts y on the other side of x and exchanges their latent parents.
            ++moved;
            ASSERT_NE(after > BOUNDARY, before > BOUNDARY);
            ASSERT_EQ(xNow, cell.yLatent);
            ASSERT_EQ(yNow, cell.xLatent);
        } else {
            ASSERT_EQ(xNow, cell.xLatent);
            ASSERT_EQ(yNow, cell.yLatent);
        }
        ++geneticsVisits[geneticsKey()];
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
    EXPECT_LT(distance(expectedGenetics, geneticsVisits), 0.005);
    EXPECT_LT(distance(expectedDurations, durationVisits), 0.005);
}

TEST(CloneOrderSwapSamplerTest, OnlyMovesAgainstClones) {
    Fixture f;
    f.x->latentGenotype(f.locus)->initializeValue(GeneticsImpl("100"));
    auto target  = std::make_shared<Target>(Target{&f});
    auto sampler = makeSampler(f, target);
    for (int i = 0; i < 100; ++i) {
        sampler.update();
    }
    EXPECT_EQ(sampler.acceptances() + sampler.rejections(), 0u);
}
