//
// Created by Maxwell Murphy on 3/10/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/containers/Infection.h"
#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/AlleleFrequenciesVector.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
#include "model/transmission_process/source_transmission_process/GeometricCOIProbability.h"

constexpr int MAX_COI = 12;
constexpr int MAX_ALLELES = 32;

TEST(MultinomialSourceTransmissionProcessTest, BasicTest) {
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using COIProbabilityImpl = GeometricCOIProbability<MAX_COI>;
    using AlleleFrequencyImpl = AlleleFrequenciesVector<MAX_ALLELES>;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequencyImpl, Locus>;
    using Infection = Infection<GeneticsImpl, Locus>;

    Locus as1("AS1", 3);
    Locus as2("AS2", 4);

    std::vector<Infection::LocusGeneticsAssignment> dlas{
            {&as1, GeneticsImpl("001")},
            {&as2, GeneticsImpl("0011")}
    };

    Infection inf1(dlas, dlas);


    std::vector<AlleleFrequencyContainer::LocusAlleleFrequencyAssignment> freqPairs{
        {&as1, AlleleFrequencyImpl({.2, .3, .5})},
        {&as2, AlleleFrequencyImpl({.1, .2, .3, .4})}
    };
    AlleleFrequencyContainer alleleFreqs(freqPairs);

    Parameter<double> coiProb(.3);
    COIProbabilityImpl coip(coiProb);

    MultinomialSourceTransmissionProcess mstp(coip, alleleFreqs, inf1);

    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb.saveState();
    EXPECT_FALSE(mstp.isDirty());
    coiProb.setValue(.95);
    EXPECT_TRUE(mstp.isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb.restoreState();
    EXPECT_FALSE(mstp.isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
}