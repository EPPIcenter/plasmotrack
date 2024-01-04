//
// Created by Maxwell Murphy on 7/18/22.
//

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/parameters/Parameter.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/SimpleLoss.h"

#include "gtest/gtest.h"


#include <fmt/core.h>

#include <memory>


using transmission_nets::core::parameters::Parameter;
using transmission_nets::core::distributions::ZTGeometric;
using transmission_nets::core::distributions::ZTPoisson;
using transmission_nets::core::containers::Locus;
using transmission_nets::core::containers::Infection;
using transmission_nets::core::containers::ParentSet;

using transmission_nets::model::transmission_process::SimpleLoss;
using transmission_nets::model::transmission_process::MultinomialSourceTransmissionProcess;

class SimpleLossTestFixture : public ::testing::Test {
    static constexpr unsigned int MAX_TRANSMISSIONS = 10;
    static constexpr unsigned int MAX_PARENTSET_SIZE = 3;
    static constexpr unsigned int MAX_ALLELES = 24;
    static constexpr unsigned int MAX_COI = 10;
    using LocusImpl = Locus;
    using GeneticsImpl = transmission_nets::core::datatypes::AllelesBitSet<MAX_ALLELES>;
    using InfectionEventImpl = Infection<GeneticsImpl, LocusImpl>;
    using AlleleFrequencyImpl = transmission_nets::core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = transmission_nets::core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl>;
    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using COIProbabilityImpl     = ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEventImpl::GenotypeParameterMap , MAX_COI>;

protected:
    SimpleLossTestFixture() : lossProb(std::make_shared<Parameter<double>>(0.9)),
                              generationProb(std::make_shared<Parameter<double>>(0.9)),
                              intp(std::make_shared<InterTransmissionProbImpl>(generationProb)),
                              simpleLoss(lossProb, intp) {

        i1->addGenetics(a1, i1_a1, i1_a1);
        i1->addGenetics(a2, i1_a2, i1_a2);
        i2->addGenetics(a1, i2_a1, i2_a1);
        i2->addGenetics(a2, i2_a2, i2_a2);
        i3->addGenetics(a1, i3_a1, i3_a1);
        i3->addGenetics(a2, i3_a2, i3_a2);
        i4->addGenetics(a1, i4_a1, i4_a1);
        i4->addGenetics(a2, i4_a2, i4_a2);

        ps.insert(i2);
        ps.insert(i3);
        ps.insert(i4);

    }

    std::shared_ptr<Parameter<double>> lossProb;
    std::shared_ptr<Parameter<double>> generationProb;
    std::shared_ptr<InterTransmissionProbImpl> intp;
    SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionImpl> simpleLoss;

    std::shared_ptr<Locus> a1 = std::make_shared<Locus>("a1", 4);
    std::shared_ptr<Locus> a2 = std::make_shared<Locus>("a2", 4);
    ParentSet<Infection<GeneticsImpl>> ps;

    std::shared_ptr<Infection<GeneticsImpl, Locus>> i1 = std::make_shared<Infection<GeneticsImpl, Locus>>("i1", 100, false);
    std::shared_ptr<Infection<GeneticsImpl, Locus>> i2 = std::make_shared<Infection<GeneticsImpl, Locus>>("i2", 100, false);
    std::shared_ptr<Infection<GeneticsImpl, Locus>> i3 = std::make_shared<Infection<GeneticsImpl, Locus>>("i3", 100, false);
    std::shared_ptr<Infection<GeneticsImpl, Locus>> i4 = std::make_shared<Infection<GeneticsImpl, Locus>>("i4", 100, false);

    GeneticsImpl i1_a1 = GeneticsImpl("111001");
    GeneticsImpl i1_a2 = GeneticsImpl("110001");

    GeneticsImpl i2_a1 = GeneticsImpl("111001");
    GeneticsImpl i2_a2 = GeneticsImpl("110001");

    GeneticsImpl i3_a1 = GeneticsImpl("111000");
    GeneticsImpl i3_a2 = GeneticsImpl("110000");

    GeneticsImpl i4_a1 = GeneticsImpl("111010");
    GeneticsImpl i4_a2 = GeneticsImpl("110010");
};

TEST_F(SimpleLossTestFixture, CoreTest) {
    long double result;
    for (size_t i = 0; i < 5000; ++i) {
        result = simpleLoss.calculateLogLikelihood(i1, ps);
    }
    fmt::print("{}\n", result);
}