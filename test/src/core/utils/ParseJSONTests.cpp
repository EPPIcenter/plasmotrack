//
// Created by Maxwell Murphy on 5/18/20.
//

#include <fstream>

#include "gtest/gtest.h"

#include "core/io/parse_json.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Alleles.h"

using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::io;

TEST(ParseJSONTests, TestMissingGenotype) {
    ASSERT_TRUE(missingGenotype("0000000"));
    ASSERT_TRUE(missingGenotype(""));
    ASSERT_FALSE(missingGenotype("00001"));
}
// A locus missing from an infection's genotype list is unobserved, like an all-zero genotype: the infection still has
// a latent genotype there, so every infection has every locus.
TEST(ParseJSONTests, UnlistedLocusIsUnobserved) {
    const auto input = nlohmann::json::parse(R"({
        "loci": [{"locus": "L0", "num_alleles": 3, "allele_freqs": [0.3, 0.3, 0.4]},
                 {"locus": "L1", "num_alleles": 2, "allele_freqs": [0.5, 0.5]}],
        "nodes": [{"id": "a", "observation_time": 400, "symptomatic": true, "allowed_parents": [],
                   "observed_genotype": [{"locus": "L0", "genotype": "010"}]}]
    })");
    auto rng        = std::make_shared<boost::random::mt19937>(1);
    const auto loci = parseLociFromJSON<Locus>(input);
    const auto infections = parseInfectionsFromJSON<Infection<AllelesBitSet<128>>, Locus>(input, 20, loci, rng);

    const auto& infection = infections.at(0);
    ASSERT_EQ(infection->loci().size(), 2u);
    EXPECT_TRUE(infection->observedGenotype().contains(loci.at("L0")));
    EXPECT_FALSE(infection->observedGenotype().contains(loci.at("L1")));
    EXPECT_EQ(infection->latentGenotype(loci.at("L1"))->value().totalPositiveCount(), 1u);
}
