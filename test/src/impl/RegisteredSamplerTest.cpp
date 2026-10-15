//
// Every sampler the Sampler set registers has an id `<type>/<target>`, unique within the set and the same for every
// chain and every Run built from the same input, so that acceptance Blocks and the checkpoint can name samplers.
// Through AbstractSampler alone, each exposes its proposal counts and its adaptation state, and both can be restored.
//

#include "gtest/gtest.h"

#include "Studies.h"

#include "impl/model/Model/Model.h"
#include "impl/model/Model/SamplerSet.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <cstdint>
#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;

namespace {

    std::vector<std::string> samplerIdsForSeed(const unsigned int seed) {
        auto rng   = std::make_shared<boost::random::mt19937>(seed);
        auto state = std::make_shared<State>(exampleStudy(), rng);
        auto model = std::make_shared<Model>(state);
        SamplerSet<Model> samplers(state, model, rng);
        return samplers.ids();
    }

    // Ids for the example input: Infections 1, 2 and 3, Loci AS1 and AS11. Every Infection has a Latent genotype at every
    // Locus, observed or not.
    const std::vector<std::string> expectedIds{
            "ConstrainedContinuousRandomWalk/mean_coi",
            "ConstrainedContinuousRandomWalk/mean_strains_tx",
            "ConstrainedContinuousRandomWalk/parent_set_size_prob",
            "ConstrainedContinuousRandomWalk/infection_duration/1",
            "RandomAllelesBitSetSampler3/latent_genotypes/1",
            "ParentRelinkSampler/infection/1",
            "StrainCountSampler/latent_genotypes/1",
            "CloneOrderSwapSampler/ordering/1",
            "InformedOrderJumpSampler/ordering/1",
            "RandomAllelesBitSetSampler4/latent_genotypes/1/AS1",
            "RandomAllelesBitSetSampler4/latent_genotypes/1/AS11",
            "ConstrainedContinuousRandomWalk/infection_duration/2",
            "RandomAllelesBitSetSampler3/latent_genotypes/2",
            "ParentRelinkSampler/infection/2",
            "StrainCountSampler/latent_genotypes/2",
            "CloneOrderSwapSampler/ordering/2",
            "InformedOrderJumpSampler/ordering/2",
            "RandomAllelesBitSetSampler4/latent_genotypes/2/AS1",
            "RandomAllelesBitSetSampler4/latent_genotypes/2/AS11",
            "ConstrainedContinuousRandomWalk/infection_duration/3",
            "RandomAllelesBitSetSampler3/latent_genotypes/3",
            "ParentRelinkSampler/infection/3",
            "StrainCountSampler/latent_genotypes/3",
            "CloneOrderSwapSampler/ordering/3",
            "InformedOrderJumpSampler/ordering/3",
            "RandomAllelesBitSetSampler4/latent_genotypes/3/AS1",
            "RandomAllelesBitSetSampler4/latent_genotypes/3/AS11",
            "AlleleSwapSampler/latent_genotypes/3/AS11",
            "ConstrainedContinuousRandomWalk/expected_false_negatives/1",
            "ConstrainedContinuousRandomWalk/expected_false_negatives/2",
            "ConstrainedContinuousRandomWalk/expected_false_negatives/3",
            "ConstrainedContinuousRandomWalk/expected_false_positives/1",
            "ConstrainedContinuousRandomWalk/expected_false_positives/2",
            "ConstrainedContinuousRandomWalk/expected_false_positives/3",
            "SALTSampler/allele_frequencies/AS1",
            "SALTSampler/allele_frequencies/AS11",
    };
}// namespace

TEST(SamplerIdsTest, TheSamplerSetNamesEverySamplerByTypeAndTarget) {
    EXPECT_EQ(samplerIdsForSeed(1), expectedIds);
}

TEST(SamplerIdsTest, IdsAreTheSameForEveryChainAndRunFromTheSameInput) {
    EXPECT_EQ(samplerIdsForSeed(1), samplerIdsForSeed(2));
}

// Analysts read ids in the acceptance Block, so a target names its Infection by the Infection's id from the input, not
// by its position.
TEST(SamplerIdsTest, TargetsNameTheInfectionByItsId) {
    auto study = exampleStudy();
    const std::map<std::string, std::string> renamed{{"1", "pf0594"}, {"2", "pf1207"}, {"3", "pf0031"}};
    for (auto& infection : study.infections) {
        infection.id = renamed.at(infection.id);
        for (auto& parent : infection.permittedParents) {
            parent = renamed.at(parent);
        }
    }
    auto rng   = std::make_shared<boost::random::mt19937>(1);
    auto state = std::make_shared<State>(study, rng);
    auto model = std::make_shared<Model>(state);
    SamplerSet<Model> samplers(state, model, rng);
    const auto& ids = samplers.ids();
    const auto registered = [&](const std::string& id) { return std::ranges::find(ids, id) != ids.end(); };

    EXPECT_TRUE(registered("ConstrainedContinuousRandomWalk/expected_false_positives/pf0594"));
    EXPECT_TRUE(registered("ConstrainedContinuousRandomWalk/expected_false_negatives/pf1207"));
    EXPECT_TRUE(registered("ConstrainedContinuousRandomWalk/infection_duration/pf0031"));
    EXPECT_TRUE(registered("RandomAllelesBitSetSampler4/latent_genotypes/pf0594/AS11"));
    // Latent parents are summed out of the model, so nothing samples them.
    EXPECT_FALSE(registered("RandomAllelesBitSetSampler/latent_parent_genotypes/pf0594_copy/AS1"));
    EXPECT_TRUE(registered("SALTSampler/allele_frequencies/AS11"));
    for (const auto& id : ids) {
        EXPECT_FALSE(id.ends_with("/1") or id.find("/1/") != std::string::npos) << id;
    }
}

// Each sampler counts the Metropolis-Hastings proposals it makes and accepts, cumulatively over the Run, so that the
// Run recorder can build acceptance Blocks as differences over each Window.
TEST(SamplerCountsTest, EachUpdateCountsOneProposalAndAnAcceptedProposalChangesTheValue) {
    auto rng   = std::make_shared<boost::random::mt19937>(7);
    auto state = std::make_shared<State>(exampleStudy(), rng);
    auto model = std::make_shared<Model>(state);
    SamplerSet<Model> samplers(state, model, rng);

    constexpr int iterations = 50;
    std::uint64_t changes    = 0;
    for (int i = 0; i < iterations; ++i) {
        const double before = state->meanCOI->value();
        samplers.step();
        changes += state->meanCOI->value() != before;
    }

    const auto meanCoi = samplers.sampler("ConstrainedContinuousRandomWalk/mean_coi").counts();
    EXPECT_EQ(meanCoi.proposed, iterations);
    EXPECT_EQ(meanCoi.accepted, changes);
    EXPECT_GT(changes, 0u);

    for (const auto& id : samplers.ids()) {
        const auto counts = samplers.sampler(id).counts();
        EXPECT_LE(counts.accepted, counts.proposed) << id;
        if (id.starts_with("SALTSampler/")) {
            // one proposal per Allele frequency visited
            EXPECT_GE(counts.proposed, iterations) << id;
        } else if (!id.starts_with("CloneOrderSwapSampler/") and !id.starts_with("InformedOrderJumpSampler/") and !id.starts_with("StrainCountSampler/")) {
            // CloneOrderSwapSampler proposes only when the Infection has a clone it can move past, and
            // InformedOrderJumpSampler only when a partner's other side has duration-prior mass, and StrainCountSampler
            // only when its remove has a second strain to remove or its add chooses a locus
            EXPECT_EQ(counts.proposed, iterations) << id;
        }
    }
}

TEST(SamplerCountsTest, RestoredCountsContinueFromWhereTheyWere) {
    auto rng   = std::make_shared<boost::random::mt19937>(7);
    auto state = std::make_shared<State>(exampleStudy(), rng);
    auto model = std::make_shared<Model>(state);
    SamplerSet<Model> samplers(state, model, rng);
    auto& meanCoi = samplers.sampler("ConstrainedContinuousRandomWalk/mean_coi");

    meanCoi.restoreCounts({.accepted = 30, .proposed = 100});
    samplers.step();

    EXPECT_EQ(meanCoi.counts().proposed, 101u);
    EXPECT_GE(meanCoi.counts().accepted, 30u);
    EXPECT_LE(meanCoi.counts().accepted, 31u);
}

namespace {
    struct Chain {
        std::shared_ptr<boost::random::mt19937> rng;
        std::shared_ptr<State> state;
        std::shared_ptr<Model> model;
        SamplerSet<Model> samplers;

        explicit Chain(const unsigned int seed, const std::vector<double>& idp = std::vector<double>(60, 1.0 / 60))
            : rng(std::make_shared<boost::random::mt19937>(seed)),
              state(std::make_shared<State>(exampleStudy(idp), rng)),
              model(std::make_shared<Model>(state)),
              samplers(state, model, rng) {}

        transmission_nets::core::samplers::AbstractSampler& sampler(const std::string& id) { return samplers.sampler(id); }
    };
}// namespace

// Samplers expose the proposal scales they have adapted to, so the final adaptation state can be recorded when burn-in
// ends. Random walks have one scale; SALT has one per Allele frequency; samplers that don't adapt have none.
TEST(AdaptationStateTest, AdaptingSamplersExposeTheirProposalScales) {
    Chain chain(11);
    for (int i = 0; i < 30; ++i) {
        chain.samplers.step();// adaptation runs from Sampling iteration 20
    }

    const auto meanCoi = chain.sampler("ConstrainedContinuousRandomWalk/mean_coi").adaptationState();
    ASSERT_EQ(meanCoi.scales.size(), 1u);
    EXPECT_EQ(meanCoi.scales[0].counts, chain.sampler("ConstrainedContinuousRandomWalk/mean_coi").counts());

    // Every adapting sampler is registered with variance 1; adaptation moves at least some of them.
    bool anyScaleMoved = false;
    for (const auto& id : chain.samplers.ids()) {
        for (const auto& scale : chain.sampler(id).adaptationState().scales) {
            anyScaleMoved = anyScaleMoved or scale.variance != 1.0;
        }
    }
    EXPECT_TRUE(anyScaleMoved);

    const auto as1 = chain.sampler("SALTSampler/allele_frequencies/AS1").adaptationState();
    EXPECT_EQ(as1.scales.size(), exampleStudy().loci[0].alleles);

    EXPECT_TRUE(chain.sampler("ParentRelinkSampler/infection/1").adaptationState().scales.empty());
}

// The checkpoint restores adaptation state into a freshly built chain; adaptation must then carry on exactly as it
// would have in the chain it was saved from.
TEST(AdaptationStateTest, RestoredStateAdaptsExactlyAsTheSavedOne) {
    Chain saved(11);
    for (int i = 0; i < 30; ++i) {
        saved.samplers.step();
    }
    Chain restored(12);
    for (const auto& id : saved.samplers.ids()) {
        restored.sampler(id).restoreAdaptationState(saved.sampler(id).adaptationState());
        EXPECT_EQ(restored.sampler(id).adaptationState(), saved.sampler(id).adaptationState()) << id;

        saved.sampler(id).adapt();
        restored.sampler(id).adapt();
        EXPECT_EQ(restored.sampler(id).adaptationState(), saved.sampler(id).adaptationState()) << id;
    }
}

TEST(AdaptationStateTest, RefusesStateWithTheWrongNumberOfScales) {
    Chain chain(11);
    const transmission_nets::core::samplers::AdaptationState twoScales{.scales = {{.variance = 1.0}, {.variance = 2.0}}};
    EXPECT_THROW(chain.sampler("ConstrainedContinuousRandomWalk/mean_coi").restoreAdaptationState(twoScales), std::invalid_argument);
    EXPECT_THROW(chain.sampler("SALTSampler/allele_frequencies/AS1").restoreAdaptationState(twoScales), std::invalid_argument);
    EXPECT_THROW(chain.sampler("ParentRelinkSampler/infection/1").restoreAdaptationState(twoScales), std::invalid_argument);
}
