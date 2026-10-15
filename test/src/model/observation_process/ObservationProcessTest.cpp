//
// Created by Maxwell Murphy on 4/19/20.
//
//
#include "gtest/gtest.h"

#include "core/computation/Accumulator.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/datatypes/Genotype.h"

#include "core/parameters/Parameter.h"

#include "model/observation_process/ObservationProcessLikelihoodv2.h"

#include <cmath>

using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;
using namespace transmission_nets::model::observation_process;

constexpr int MAX_ALLELES = 32;

TEST(ObservationProcessTest, V2PairsEachCountWithItsRate) {

    // Latent 110000 observed as 101000: one true positive, one false negative, one false positive, three true negatives.
    auto observed = std::make_shared<Data<Genotype>>(Genotype("101000"));
    auto latent   = std::make_shared<Parameter<Genotype>>(Genotype("110000"));

    // Expected errors per locus, scaled by the 6 alleles into per-allele rates of .05 and .2.
    auto expectedFalsePositives = std::make_shared<Parameter<double>>(.3);
    auto expectedFalseNegatives = std::make_shared<Parameter<double>>(1.2);

    ObservationProcessLikelihoodv2 target(observed, latent, expectedFalsePositives, expectedFalseNegatives);

    const double fpr = .05;
    const double fnr = .2;
    EXPECT_DOUBLE_EQ(target.value(), std::log(1 - fnr) + 3 * std::log(1 - fpr) + std::log(fpr) + std::log(fnr));
}
