//
// The Study and the State built from it: the JSON adapter reads every field and nothing else; the State applies the
// model's rules to the data (an unobserved Locus gets one random Allele, a starting genotype is capped at the COI
// limit, frequencies are floored, the null model ignores genotypes and times) and refuses a Study that doesn't fit
// itself.
//

#include "gtest/gtest.h"

#include "Studies.h"

#include "impl/model/Model/Refusal.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;

namespace {
    using Engine = boost::random::mt19937;

    // One Locus of three Alleles; "a" observed 010 there, "b" observed nothing, "b" may be infected by "a".
    Study twoInfections() {
        return {.loci                            = {uniformLocus("L0", 3)},
                .infections                      = {{.id = "a", .observationTime = 400, .observedGenotypes = {{"L0", "010"}}},
                                                    {.id = "b", .observationTime = 420, .symptomatic = false, .permittedParents = {"a"}}},
                .symptomaticInfectionDurations   = flatDurations(),
                .asymptomaticInfectionDurations  = flatDurations()};
    }

    std::shared_ptr<State> stateOf(const Study& study, const bool nullModel = false) {
        return std::make_shared<State>(study, std::make_shared<Engine>(1), nullModel);
    }
}// namespace

TEST(StudyTest, TheJsonAdapterReadsEveryFieldAndDefaultsSymptomatic) {
    const auto input = nlohmann::json::parse(R"({
        "loci": [{"locus": "L0", "num_alleles": 3, "allele_freqs": [0.3, 0.3, 0.4]},
                 {"locus": "L1", "num_alleles": 2, "allele_freqs": [0.5, 0.5]}],
        "nodes": [{"id": "a", "observation_time": 400, "symptomatic": false, "allowed_parents": [],
                   "observed_genotype": [{"locus": "L1", "genotype": "01"}, {"locus": "L0", "genotype": "010"}]},
                  {"id": "b", "observation_time": 420.5, "allowed_parents": ["a"], "observed_genotype": []}]
    })");

    const auto study = studyFromJson(input, flatDurations(30), flatDurations(45));

    ASSERT_EQ(study.loci.size(), 2u);
    EXPECT_EQ(study.loci[0].label, "L0");
    EXPECT_EQ(study.loci[0].alleles, 3u);
    EXPECT_EQ(study.loci[1].alleleFrequencies, (std::vector<double>{0.5, 0.5}));
    ASSERT_EQ(study.infections.size(), 2u);
    EXPECT_EQ(study.infections[0].id, "a");
    EXPECT_EQ(study.infections[0].observationTime, 400);
    EXPECT_FALSE(study.infections[0].symptomatic);
    ASSERT_EQ(study.infections[0].observedGenotypes.size(), 2u);
    EXPECT_EQ(study.infections[0].observedGenotypes[0].locus, "L1");// input order, not Locus order
    EXPECT_EQ(study.infections[0].observedGenotypes[1].alleles, "010");
    EXPECT_TRUE(study.infections[1].symptomatic);
    EXPECT_EQ(study.infections[1].observationTime, 420.5);
    EXPECT_EQ(study.infections[1].permittedParents, std::vector<std::string>{"a"});
    EXPECT_EQ(study.symptomaticInfectionDurations.size(), 30u);
    EXPECT_EQ(study.asymptomaticInfectionDurations.size(), 45u);
}

TEST(StudyTest, TheJsonAdapterRefusesADocumentThatLacksAField) {
    const auto noNodes = nlohmann::json::parse(R"({"loci": []})");
    EXPECT_THROW(studyFromJson(noNodes, {}, {}), Refused);
    const auto noTime = nlohmann::json::parse(R"({"loci": [], "nodes": [{"id": "a", "allowed_parents": [], "observed_genotype": []}]})");
    EXPECT_THROW(studyFromJson(noTime, {}, {}), Refused);
}

// An Infection has a Latent genotype at every Locus of the Study. Where nothing was observed, whether the Locus is
// unlisted or its genotype all zeros, it starts as one random Allele; where something was, as the observation.
TEST(StudyTest, AnUnobservedLocusStartsAsOneRandomAllele) {
    auto study = twoInfections();
    study.loci.push_back(uniformLocus("L1", 4));
    study.infections[0].observedGenotypes.push_back({"L1", "0000"});
    const auto state = stateOf(study);

    const auto a  = state->infection("a");
    const auto b  = state->infection("b");
    const auto l0 = state->loci.at("L0");
    const auto l1 = state->loci.at("L1");
    EXPECT_EQ(a->loci().size(), 2u);
    EXPECT_EQ(b->loci().size(), 2u);
    EXPECT_TRUE(a->observedGenotype().contains(l0));
    EXPECT_FALSE(a->observedGenotype().contains(l1));
    EXPECT_FALSE(b->observedGenotype().contains(l0));
    EXPECT_EQ(a->latentGenotype(l0)->value().allelesStr(), "010");
    EXPECT_EQ(a->latentGenotype(l1)->value().totalPositiveCount(), 1u);
    EXPECT_EQ(b->latentGenotype(l0)->value().totalPositiveCount(), 1u);
    EXPECT_EQ(b->latentGenotype(l1)->value().totalPositiveCount(), 1u);
}

// A starting genotype holds fewer Alleles than the COI limit: the lowest-numbered of those observed.
TEST(StudyTest, AStartingGenotypeIsCappedAtTheCOILimit) {
    auto study = twoInfections();
    study.loci[0] = uniformLocus("L0", 30);
    study.infections[0].observedGenotypes = {{"L0", std::string(30, '1')}};
    const auto state = stateOf(study);

    const auto latent = state->infection("a")->latentGenotype(state->loci.at("L0"))->value();
    EXPECT_LT(latent.totalPositiveCount(), static_cast<unsigned int>(MAX_COI));
    EXPECT_GT(latent.totalPositiveCount(), 0u);
    for (unsigned int allele = 0; allele < 30; ++allele) {
        EXPECT_EQ(latent.allele(allele), allele < latent.totalPositiveCount()) << "allele " << allele;
    }
}

// Frequencies below the minimum are raised to it before the Simplex renormalises, so none starts at zero or far below
// the minimum.
TEST(StudyTest, AlleleFrequenciesAreFlooredAtTheMinimum) {
    auto study = twoInfections();
    study.loci[0].alleleFrequencies = {0.0, 1e-9, 1.0};
    const auto state = stateOf(study);

    const auto frequencies = state->alleleFrequencies->alleleFrequencies(state->loci.at("L0"))->value().frequencies();
    ASSERT_EQ(frequencies.size(), 3u);
    for (const auto frequency : frequencies) {
        EXPECT_NEAR(frequency, std::max(frequency, MIN_ALLELE_FREQUENCY), 1e-5 * MIN_ALLELE_FREQUENCY);
    }
    EXPECT_GT(frequencies[0], 0.0);
}

TEST(StudyTest, TheNullModelIgnoresGenotypesAndObservationTimes) {
    const auto state = stateOf(twoInfections(), true);

    for (const auto& infection : state->infections) {
        EXPECT_TRUE(infection->observedGenotype().empty()) << infection->id();
        EXPECT_EQ(infection->observationTime()->value(), 1000) << infection->id();
        EXPECT_EQ(infection->latentGenotype(state->loci.at("L0"))->value().totalPositiveCount(), 1u) << infection->id();
    }
}

TEST(StudyTest, InfectionsAreFoundByIdAndByPosition) {
    const auto state = stateOf(twoInfections());

    EXPECT_EQ(state->indexOf("b"), 1u);
    EXPECT_EQ(state->infection("b"), state->infections[1]);
    EXPECT_THROW(static_cast<void>(state->indexOf("c")), std::out_of_range);
    EXPECT_THROW(static_cast<void>(state->infection("c")), std::out_of_range);
}

TEST(StudyTest, AStudyThatDoesNotFitItselfIsRefused) {
    auto unknownLocus = twoInfections();
    unknownLocus.infections[0].observedGenotypes = {{"L9", "010"}};
    EXPECT_THROW(stateOf(unknownLocus), Refused);

    auto unknownParent = twoInfections();
    unknownParent.infections[1].permittedParents = {"z"};
    EXPECT_THROW(stateOf(unknownParent), Refused);

    auto wrongGenotype = twoInfections();
    wrongGenotype.infections[0].observedGenotypes = {{"L0", "0100"}};
    EXPECT_THROW(stateOf(wrongGenotype), Refused);

    auto wrongFrequencies = twoInfections();
    wrongFrequencies.loci[0].alleleFrequencies = {0.5, 0.5};
    EXPECT_THROW(stateOf(wrongFrequencies), Refused);

    auto negativeTime = twoInfections();
    negativeTime.infections[0].observationTime = -1;
    EXPECT_THROW(stateOf(negativeTime), Refused);

    auto duplicateLocus = twoInfections();
    duplicateLocus.loci.push_back(uniformLocus("L0", 3));
    EXPECT_THROW(stateOf(duplicateLocus), Refused);
}
