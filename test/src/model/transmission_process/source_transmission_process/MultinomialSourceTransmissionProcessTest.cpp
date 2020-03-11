//
// Created by Maxwell Murphy on 3/10/20.
//

#include "gtest/gtest.h"

#include "core/containers/Infection.h"
#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/AlleleFrequenciesVector.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
#include "model/transmission_process/source_transmission_process/GeometricCOIProbability.h"

constexpr int MAX_COI = 1200;
constexpr int MAX_ALLELES = 32;

TEST(MultinomialSourceTransmissionProcessTest, BasicTest) {
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using COIProbabilityImpl = GeometricCOIProbability<MAX_COI>;
    using AlleleFrequencyImpl = AlleleFrequenciesVector<MAX_ALLELES>;

    Locus as1("AS1", 3);
    Locus as2("AS2", 4);

    std::vector<LocusGeneticsAssignment<GeneticsImpl, Locus>> dlas{
            {&as1, GeneticsImpl("001")},
            {&as2, GeneticsImpl("0011")}
    };

    Infection<GeneticsImpl> inf1(dlas, dlas);


    std::vector<LocusAlleleFrequencyAssignment<AlleleFrequencyImpl, Locus>> freqPairs{
        {&as1, AlleleFrequencyImpl({.2, .3, .5})},
        {&as2, AlleleFrequencyImpl({.1, .2, .3, .4})}
    };
    AlleleFrequencyContainer<AlleleFrequencyImpl> alleleFreqs(freqPairs);

    Parameter<double> coiProb(.3);
    COIProbabilityImpl coip(coiProb);

    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyImpl> mstp(coip, alleleFreqs);

    std::cout << "LogLikelihood: " << mstp.calculateLikelihood(inf1) << std::endl;
    coiProb.saveState();
    coiProb.setValue(.95);
    std::cout << "LogLikelihood: " << mstp.calculateLikelihood(inf1) << std::endl;
    coiProb.restoreState();
    std::cout << "LogLikelihood: " << mstp.calculateLikelihood(inf1) << std::endl;
}