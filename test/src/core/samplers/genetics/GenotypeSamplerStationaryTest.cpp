//
// A Metropolis-Hastings sampler must leave its target invariant: the states it visits are visited in proportion to
// the target. Each sampler here runs alone against a target that weights each allele independently, and the visited
// genotypes' frequencies are compared with the target restricted to the genotypes visited.
//

#include "gtest/gtest.h"

#include "core/samplers/genetics/AlleleSwapSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/specialized/ParentRelinkSampler.h"
#include "impl/model/Model/config.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::parameters::Parameter;

// The parent-set samplers below propose with error rates of 0.3 rather than their defaults: with rates near 0.5
// every genotype's proposal probability is large, so frequencies settle quickly, whereas rarely proposed genotypes
// are also left rarely and their frequencies converge slowly.

namespace {
    using Engine = boost::random::mt19937;

    constexpr std::array<double, 3> WEIGHTS{0.7, -0.4, 0.2};

    // log target = sum of the weights of the alleles present.
    struct WeightedAllelesTarget {
        std::shared_ptr<Parameter<GeneticsImpl>> genotype;

        double value() const {
            double out = 0;
            for (unsigned int j = 0; j < WEIGHTS.size(); ++j) {
                out += genotype->value().allele(j) ? WEIGHTS[j] : 0.0;
            }
            return out;
        }

        [[nodiscard]] bool isDirty() const { return false; }
    };

    double targetWeight(const std::string& genotype) {
        return std::exp(WeightedAllelesTarget{std::make_shared<Parameter<GeneticsImpl>>(GeneticsImpl(genotype))}.value());
    }

    // Total variation distance between the visited genotypes' frequencies and the target over those genotypes, each
    // weighted by `support`, the measure of the rest of the state over which the sampler can reach it.
    template<typename Support>
    double distanceFromTarget(const std::map<std::string, int>& visits, const int draws, Support&& support) {
        double total = 0;
        for (const auto& [genotype, count] : visits) {
            total += targetWeight(genotype) * support(genotype);
        }
        double distance = 0;
        for (const auto& [genotype, count] : visits) {
            const double expected = targetWeight(genotype) * support(genotype) / total;
            distance += std::abs(static_cast<double>(count) / draws - expected) / 2;
        }
        return distance;
    }

    double distanceFromTarget(const std::map<std::string, int>& visits, const int draws) {
        return distanceFromTarget(visits, draws, [](const std::string&) { return 1.0; });
    }

    template<typename Sampler>
    std::map<std::string, int> run(Sampler& sampler, const std::shared_ptr<Parameter<GeneticsImpl>>& genotype, const int draws) {
        for (int i = 0; i < 5000; ++i) {
            sampler.update();
        }
        std::map<std::string, int> visits;
        for (int i = 0; i < draws; ++i) {
            sampler.update();
            ++visits[genotype->value().serialize()];
        }
        return visits;
    }

    // One three-allele locus. "c" may be infected by "p"; "q" has no candidate parents; "r" may be infected by "p"
    // only while its infection duration is under 30 days, since p was infected on day 390; "m" may be infected by "p"
    // and has no observation at the locus.
    nlohmann::json makeInput() {
        nlohmann::json input;
        input["loci"].push_back({{"locus", "L0"}, {"num_alleles", 3}, {"allele_freqs", {1.0 / 3, 1.0 / 3, 1.0 / 3}}});
        auto node = [&](const std::string& id, const double time, const std::string& observed, const std::vector<std::string>& parents) {
            input["nodes"].push_back({{"id", id},
                                      {"observation_time", time},
                                      {"symptomatic", true},
                                      {"observed_genotype", {{{"locus", "L0"}, {"genotype", observed}}}},
                                      {"allowed_parents", parents}});
        };
        node("p", 400, "110", {});
        node("c", 500, "100", {"p"});
        node("q", 450, "010", {});
        node("r", 420, "010", {"p"});
        node("m", 500, "000", {"p"});// no observation at the locus
        return input;
    }

    struct Fixture {
        std::shared_ptr<Engine> rng = std::make_shared<Engine>(20261008);
        std::shared_ptr<State> state;
        std::shared_ptr<LocusImpl> locus;

        Fixture() {
            const std::vector<double> idp(60, 1.0 / 60);
            state = std::make_shared<State>(makeInput(), idp, idp, rng);
            locus = state->loci.at("L0");
        }

        size_t index(const std::string& id) const {
            for (size_t i = 0; i < state->infections.size(); ++i) {
                if (state->infections[i]->id() == id) {
                    return i;
                }
            }
            throw std::out_of_range(id);
        }

        void setGenotype(const std::shared_ptr<InfectionEvent>& infection, const std::string& genotype) const {
            infection->latentGenotype(locus)->initializeValue(GeneticsImpl(genotype));
        }
    };
}// namespace

TEST(GenotypeSamplerStationaryTest, ParentSetGenotypeSampler) {
    Fixture f;
    const auto c = f.index("c");
    f.setGenotype(f.state->infections[f.index("p")], "110");
    f.setGenotype(f.state->latentParents[c], "011");
    f.setGenotype(f.state->infections[c], "100");

    const auto& child = f.state->infections[c];
    auto parentSet    = std::make_shared<ParentSetImpl>(f.state->infectionEventOrdering, child, f.state->allowedRelationships->allowedParents(child));
    ASSERT_EQ(parentSet->value().size(), 1u);

    auto genotype = child->latentGenotype(f.locus);
    auto target   = std::make_shared<WeightedAllelesTarget>(WeightedAllelesTarget{genotype});
    transmission_nets::core::samplers::genetics::RandomAllelesBitSetSampler3<WeightedAllelesTarget, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>
            sampler(child, parentSet, f.state->latentParents[c], target, f.rng, 0.3, 0.3);

    constexpr int draws = 400000;
    const auto visits   = run(sampler, genotype, draws);
    EXPECT_GE(visits.size(), 3u);
    EXPECT_LT(distanceFromTarget(visits, draws), 0.01);
}

TEST(GenotypeSamplerStationaryTest, ParentSetGenotypeSamplerAtUnobservedLocus) {
    Fixture f;
    const auto m = f.index("m");
    f.setGenotype(f.state->infections[f.index("p")], "110");
    f.setGenotype(f.state->latentParents[m], "011");
    f.setGenotype(f.state->infections[m], "100");

    const auto& child = f.state->infections[m];
    ASSERT_FALSE(child->observedGenotype().contains(f.locus));
    auto parentSet = std::make_shared<ParentSetImpl>(f.state->infectionEventOrdering, child, f.state->allowedRelationships->allowedParents(child));
    ASSERT_EQ(parentSet->value().size(), 1u);

    auto genotype = child->latentGenotype(f.locus);
    auto target   = std::make_shared<WeightedAllelesTarget>(WeightedAllelesTarget{genotype});
    transmission_nets::core::samplers::genetics::RandomAllelesBitSetSampler3<WeightedAllelesTarget, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>
            sampler(child, parentSet, f.state->latentParents[m], target, f.rng, 0.3, 0.3);

    constexpr int draws = 400000;
    const auto visits   = run(sampler, genotype, draws);
    EXPECT_EQ(visits.size(), 7u);// every genotype with at least one allele
    EXPECT_LT(distanceFromTarget(visits, draws), 0.01);
}

// One try is a plain Metropolis-Hastings independence move; more tries make it a multiple-try move.
class ParentRelinkStationaryTest : public testing::TestWithParam<unsigned int> {};

INSTANTIATE_TEST_SUITE_P(Tries, ParentRelinkStationaryTest, testing::Values(1u, 4u), [](const auto& info) { return std::to_string(info.param); });

TEST_P(ParentRelinkStationaryTest, LeavesTargetInvariant) {
    Fixture f;
    const auto r = f.index("r");
    f.setGenotype(f.state->infections[f.index("p")], "110");
    f.setGenotype(f.state->latentParents[r], "011");
    f.setGenotype(f.state->infections[r], "010");

    const auto& infection = f.state->infections[r];
    auto genotype         = infection->latentGenotype(f.locus);
    auto target           = std::make_shared<WeightedAllelesTarget>(WeightedAllelesTarget{genotype});
    // Durations are proposed from a shape unlike the (flat) target, which the correction must undo.
    std::vector<double> shape(60);
    for (size_t day = 0; day < shape.size(); ++day) {
        shape[day] = 1.0 + static_cast<double>(day);
    }
    transmission_nets::core::samplers::specialized::ParentRelinkSampler<WeightedAllelesTarget, Engine, InfectionEvent>
            sampler(infection, f.state->latentParents[r], f.state->allowedRelationships->allowedParents(infection), target, f.rng, 1.0, 60.0, shape, 0.3, 0.3, GetParam());

    // A genotype can be drawn from a source if it lies within the source's alleles and shares one of them, at the
    // durations for which the source came first: [1, 30) for p, infected on day 390, and for the latent parent any
    // duration whose day, round(duration), is below 60, so [1, 59.5).
    const std::vector<std::pair<GeneticsImpl, double>> sources{{GeneticsImpl("110"), 29.0}, {GeneticsImpl("011"), 58.5}};
    auto durations = [&](const std::string& genotypeString) {
        const GeneticsImpl candidate(genotypeString);
        double longest = 0;
        for (const auto& [alleles, length] : sources) {
            if (GeneticsImpl::falsePositiveCount(alleles, candidate) == 0 and GeneticsImpl::truePositiveCount(alleles, candidate) > 0) {
                longest = std::max(longest, length);
            }
        }
        return longest;
    };

    for (int i = 0; i < 5000; ++i) {
        sampler.update();
    }
    std::map<std::string, int> visits;
    constexpr int draws = 400000;
    int early           = 0;
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        ++visits[genotype->value().serialize()];
        early += infection->infectionDuration()->value() < 30.0 ? 1 : 0;
    }
    double total      = 0;
    double totalEarly = 0;
    for (const auto& [genotypeString, count] : visits) {
        total += targetWeight(genotypeString) * durations(genotypeString);
        totalEarly += targetWeight(genotypeString) * 29.0;
    }
    EXPECT_EQ(visits.size(), 5u);// every genotype within one source's alleles
    EXPECT_NEAR(static_cast<double>(early) / draws, totalEarly / total, 0.01);
    EXPECT_LT(distanceFromTarget(visits, draws, durations), 0.01);
}

// The swap keeps the number of alleles, so from a one-allele genotype it visits exactly the one-allele genotypes, in
// proportion to the target.
TEST(GenotypeSamplerStationaryTest, AlleleSwapSampler) {
    Fixture f;
    const auto c = f.index("c");
    f.setGenotype(f.state->infections[c], "100");
    auto genotype = f.state->infections[c]->latentGenotype(f.locus);
    auto target   = std::make_shared<WeightedAllelesTarget>(WeightedAllelesTarget{genotype});
    transmission_nets::core::samplers::genetics::AlleleSwapSampler<WeightedAllelesTarget, Engine, GeneticsImpl> sampler(genotype, target, f.rng);

    constexpr int draws = 200000;
    const auto visits   = run(sampler, genotype, draws);
    EXPECT_EQ(visits.size(), 3u);
    EXPECT_LT(distanceFromTarget(visits, draws), 0.01);
}
