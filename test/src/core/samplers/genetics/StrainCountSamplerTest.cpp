//
// StrainCountSampler must leave its target invariant. One infection has two three-allele loci, one with an observed
// allele and one with nothing observed; the target is an arbitrary function of both genotypes. The move alone does not
// reach every genotype pair, so rather than running it as a chain, exact draws from the target are each moved once and
// must still follow the target.
//

#include "gtest/gtest.h"

#include "core/samplers/genetics/StrainCountSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"
#include "../../../impl/Studies.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;
using transmission_nets::core::samplers::genetics::StrainCountSampler;

namespace {
    using Engine = boost::random::mt19937;

    const std::vector<std::string> GENOTYPES{"100", "010", "001", "110", "101", "011", "111"};

    // An arbitrary weight for each locus' genotype, and a term coupling the loci by their allele counts.
    double logTarget(const Genotype& a, const Genotype& b) {
        const std::map<std::string, double> wa{{"100", 0.3}, {"010", -0.2}, {"001", 0.5}, {"110", 0.1}, {"101", -0.4}, {"011", 0.7}, {"111", -0.1}};
        const std::map<std::string, double> wb{{"100", -0.5}, {"010", 0.4}, {"001", 0.0}, {"110", 0.6}, {"101", 0.2}, {"011", -0.3}, {"111", 0.4}};
        const int ca = static_cast<int>(a.totalPositiveCount()), cb = static_cast<int>(b.totalPositiveCount());
        return wa.at(a.serialize()) + wb.at(b.serialize()) - 0.4 * std::abs(ca - cb);
    }

    // One Infection with two Loci: two Alleles observed at L0, nothing at L1.
    Study makeStudy() {
        return {.loci                           = {uniformLocus("L0", 3), uniformLocus("L1", 3)},
                .infections                     = {{.id = "x", .observationTime = 400, .observedGenotypes = {{"L0", "011"}, {"L1", "000"}}}},
                .symptomaticInfectionDurations  = flatDurations(),
                .asymptomaticInfectionDurations = flatDurations()};
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261012);
        std::shared_ptr<State> state;
        std::shared_ptr<InfectionEvent> x;
        std::shared_ptr<LocusImpl> l0, l1;

        Fixture() {
            state = std::make_shared<State>(makeStudy(), rng);
            x     = state->infections[0];
            l0    = state->loci.at("L0");
            l1    = state->loci.at("L1");
        }
    };

    struct Target {
        const Fixture* f;
        double value() const { return logTarget(f->x->latentGenotype(f->l0)->value(), f->x->latentGenotype(f->l1)->value()); }
        [[nodiscard]] bool isDirty() const { return false; }
    };
}// namespace

class StrainCountInvarianceTest : public testing::TestWithParam<unsigned int> {};

// The largest number of alleles a locus may hold: 3 allows every genotype; 2 makes adds to two alleles impossible.
INSTANTIATE_TEST_SUITE_P(MaxCOI, StrainCountInvarianceTest, testing::Values(3u, 2u), [](const auto& info) { return std::to_string(info.param); });

TEST_P(StrainCountInvarianceTest, LeavesTargetInvariant) {
    Fixture f;
    auto target = std::make_shared<Target>(Target{&f});
    StrainCountSampler<Target, Engine, InfectionEvent> sampler(f.x, target, f.rng, GetParam(), 0.5, 0.2, 4.0);

    // Exact draws over genotype pairs allowed by maxCOI.
    std::vector<std::pair<std::string, std::string>> states;
    std::vector<double> cumulative;
    double total = 0;
    for (const auto& a : GENOTYPES) {
        for (const auto& b : GENOTYPES) {
            if (Genotype(a).totalPositiveCount() > GetParam() or Genotype(b).totalPositiveCount() > GetParam()) {
                continue;
            }
            total += std::exp(logTarget(Genotype(a), Genotype(b)));
            states.emplace_back(a, b);
            cumulative.push_back(total);
        }
    }

    boost::random::uniform_01<> uniform;
    constexpr int draws = 400000;
    std::map<std::pair<std::string, std::string>, int> visits;
    int moved = 0;
    for (int i = 0; i < draws; ++i) {
        const double u     = uniform(*f.rng) * total;
        const auto& [a, b] = states[std::min(states.size() - 1, static_cast<size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()))];
        f.x->latentGenotype(f.l0)->initializeValue(Genotype(a));
        f.x->latentGenotype(f.l1)->initializeValue(Genotype(b));
        sampler.update();
        const auto after = std::pair{f.x->latentGenotype(f.l0)->value().serialize(), f.x->latentGenotype(f.l1)->value().serialize()};
        moved += after != std::pair{a, b};
        ASSERT_LE(Genotype(after.first).totalPositiveCount(), GetParam());
        ASSERT_LE(Genotype(after.second).totalPositiveCount(), GetParam());
        ++visits[after];
    }

    double distance = 0;
    for (size_t i = 0; i < states.size(); ++i) {
        const double p = std::exp(logTarget(Genotype(states[i].first), Genotype(states[i].second))) / total;
        distance += std::abs((visits.contains(states[i]) ? visits.at(states[i]) : 0) / static_cast<double>(draws) - p) / 2;
    }
    EXPECT_GT(moved, draws / 20);
    EXPECT_LT(distance, 0.01);
}
