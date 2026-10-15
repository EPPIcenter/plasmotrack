//
// Under a flat target every latent genotype is equally likely, so a correct Metropolis-Hastings sampler must
// leave each allele present half of the time, whatever its proposal favours.
//

#include "gtest/gtest.h"

#include "core/containers/AllowedRelationships.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Genotype.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler4.h"

#include <boost/random.hpp>

#include <array>
#include <memory>
#include <string>

using namespace transmission_nets::core;

namespace {
    using Genotype  = datatypes::Genotype;
    using Infection = containers::Infection<>;

    struct FlatTarget {
        double value() { return 0.0; }
    };
}// namespace

TEST(RandomAllelesBitSetSampler4Test, FlatTargetLeavesEveryAlleleEquallyLikely) {
    auto locus = std::make_shared<containers::Locus>("L0", 2);

    // A path a -> b -> c, so proposals walk through children and each step's proposal depends on the step before.
    auto a = std::make_shared<Infection>("a", 10.0);
    auto b = std::make_shared<Infection>("b", 20.0);
    auto c = std::make_shared<Infection>("c", 30.0);
    a->addGenetics(locus, Genotype("10"), Genotype("10"));
    b->addGenetics(locus, Genotype("10"), Genotype("10"));
    c->addGenetics(locus, Genotype("01"), Genotype("01"));
    auto latentParent = std::make_shared<Infection>(*a);

    auto relationships = std::make_shared<containers::AllowedRelationships<Infection>>();
    for (const auto& infection : {a, b, c, latentParent}) {
        relationships->addInfectionEvent(infection);
    }
    relationships->addChild(a, b);
    relationships->addChild(b, c);

    auto rng    = std::make_shared<boost::random::mt19937>(20261008);
    auto target = std::make_shared<FlatTarget>();
    samplers::genetics::RandomAllelesBitSetSampler4<FlatTarget, boost::random::mt19937, Infection, containers::Locus>
            sampler(a, latentParent, locus, relationships, target, rng, 20);

    for (int i = 0; i < 10000; ++i) {
        sampler.update();
    }

    const std::array<std::shared_ptr<Infection>, 3> path{a, b, c};
    std::array<std::array<double, 2>, 3> present{};
    constexpr int draws = 400000;
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (size_t n = 0; n < path.size(); ++n) {
            for (unsigned int j = 0; j < 2; ++j) {
                present[n][j] += path[n]->latentGenotype(locus)->value().allele(j);
            }
        }
    }

    for (size_t n = 0; n < path.size(); ++n) {
        for (unsigned int j = 0; j < 2; ++j) {
            EXPECT_NEAR(present[n][j] / draws, 0.5, 0.02) << "infection " << path[n]->id() << ", allele " << j;
        }
    }
}
