//
// Created by Maxwell Murphy on 3/10/20.
//
//
#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/distributions/ZTGeometric.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

constexpr int MAX_COI     = 12;
constexpr int MAX_ALLELES = 32;

TEST(MultinomialSourceTransmissionProcessTest, BasicTest) {
    using GeneticsImpl             = AllelesBitSet<MAX_ALLELES>;
    using COIProbabilityImpl       = ZTGeometric<MAX_COI>;
    using AlleleFrequencyImpl      = Simplex;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequencyImpl, Locus>;
    using Infection                = Infection<GeneticsImpl, Locus>;

    auto as1 = std::make_shared<Locus>("AS1", 3);
    auto as2 = std::make_shared<Locus>("AS2", 4);

    auto inf1 = std::make_shared<Infection>("inf1", 10.0);
    inf1->addGenetics(as1, "001", "001");
    inf1->addGenetics(as2, "0011", "0011");

    auto alleleFreqs = std::make_shared<AlleleFrequencyContainer>();
    alleleFreqs->addLocus(as1);
    alleleFreqs->alleleFrequencies(as1)->initializeValue({.2, .3, .5});
    alleleFreqs->addLocus(as2);
    alleleFreqs->alleleFrequencies(as2)->initializeValue({.1, .2, .3, .4});

    auto coiProb = std::make_shared<Parameter<double>>(.3);
    auto coip    = std::make_shared<COIProbabilityImpl>(coiProb);

    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, Infection::GenotypeParameterMap, MAX_COI> mstp(coip, alleleFreqs, inf1->loci(), inf1->latentGenotype());

    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb->saveState("state1");
    EXPECT_FALSE(mstp.isDirty());
    coiProb->setValue(.001);
    EXPECT_TRUE(mstp.isDirty());
    EXPECT_TRUE(coip->isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb->restoreState("state1");
    EXPECT_FALSE(mstp.isDirty());
    EXPECT_FALSE(coip->isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
}