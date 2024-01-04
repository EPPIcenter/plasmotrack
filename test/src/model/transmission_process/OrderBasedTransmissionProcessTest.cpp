//
// Created by Maxwell Murphy on 3/6/20.
//
//
#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Ordering.h"

#include "core/computation/OrderDerivedParentSet.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionNoMutation.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

constexpr int MAX_PARENTS       = 1;
constexpr int MAX_ALLELES       = 32;
constexpr int MAX_COI           = 10;
constexpr int MAX_TRANSMISSIONS = 5;

TEST(OrderBasedTransmissionProcessTest, CoreTest) {
    using GeneticsImpl             = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent           = Infection<GeneticsImpl>;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<Simplex>;
    using Ordering                 = Ordering<InfectionEvent>;

    using COITransitionProbImpl     = ZTMultiplicativeBinomial<MAX_COI>;
    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl      = NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>;

    using COIProbabilityImpl     = ZTGeometric<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEvent::GenotypeParameterMap, MAX_COI>;
    using TransmissionProcess    = OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent, OrderDerivedParentSet<InfectionEvent, Ordering>>;

    auto as1 = std::make_shared<Locus>("AS1", 5);
    auto as2 = std::make_shared<Locus>("AS2", 6);

    auto inf1 = std::make_shared<InfectionEvent>("1", 100, false);
    auto inf2 = std::make_shared<InfectionEvent>("2", 100, false);
    auto inf3 = std::make_shared<InfectionEvent>("3", 100, false);
    auto inf4 = std::make_shared<InfectionEvent>("4", 100, false);

    inf1->addGenetics(as1, "11010", "11010");
    inf1->addGenetics(as2, "000011", "000011");
    inf2->addGenetics(as1, "11010", "11010");
    inf2->addGenetics(as2, "000011", "000011");
    inf3->addGenetics(as1, "11010", "11010");
    inf3->addGenetics(as2, "000011", "000011");
    inf4->addGenetics(as1, "11010", "11010");
    inf4->addGenetics(as2, "000011", "000011");

    auto infectionOrder = std::make_shared<Ordering>(std::vector{inf1, inf2, inf3, inf4});

    auto ps1 = std::make_shared<OrderDerivedParentSet<InfectionEvent, Ordering>>(infectionOrder, inf1, std::vector{inf2, inf3, inf4});
    auto ps2 = std::make_shared<OrderDerivedParentSet<InfectionEvent, Ordering>>(infectionOrder, inf2, std::vector{inf1, inf3, inf4});
    auto ps3 = std::make_shared<OrderDerivedParentSet<InfectionEvent, Ordering>>(infectionOrder, inf3, std::vector{inf1, inf2, inf4});
    auto ps4 = std::make_shared<OrderDerivedParentSet<InfectionEvent, Ordering>>(infectionOrder, inf4, std::vector{inf1, inf2, inf3});

    //    std::vector<AlleleFrequencyContainer::LocusAlleleFrequencyAssignment> lfas {
    //            {&as1, Simplex(as1.totalAlleles())},
    //            {&as2, Simplex(as2.totalAlleles())}
    //    };

    auto afc = std::make_shared<AlleleFrequencyContainer>();
    afc->addLocus(as1);
    afc->addLocus(as2);

    auto ztmbProb  = std::make_shared<Parameter<double>>(.3);
    auto ztmbAssoc = std::make_shared<Parameter<double>>(1.0);
    auto ztmbCTP   = std::make_shared<COITransitionProbImpl>(ztmbProb, ztmbAssoc);

    auto geoGenProb = std::make_shared<Parameter<double>>(.8);
    auto geoGen     = std::make_shared<InterTransmissionProbImpl>(geoGenProb);

    auto nodeTransmission = std::make_shared<NodeTransmissionImpl>(ztmbCTP, geoGen);

    auto geoCOIProb = std::make_shared<Parameter<double>>(.9);
    auto geoCOI     = std::make_shared<COIProbabilityImpl>(geoCOIProb);

    auto mstp1 = std::make_shared<SourceTransmissionImpl>(geoCOI, afc, inf1->loci(), inf1->latentGenotype());
    auto mstp2 = std::make_shared<SourceTransmissionImpl>(geoCOI, afc, inf2->loci(), inf2->latentGenotype());
    auto mstp3 = std::make_shared<SourceTransmissionImpl>(geoCOI, afc, inf3->loci(), inf3->latentGenotype());
    auto mstp4 = std::make_shared<SourceTransmissionImpl>(geoCOI, afc, inf4->loci(), inf4->latentGenotype());

    auto tp1 = std::make_shared<TransmissionProcess>(nodeTransmission, mstp1, inf1, ps1);
    auto tp2 = std::make_shared<TransmissionProcess>(nodeTransmission, mstp2, inf2, ps2);
    auto tp3 = std::make_shared<TransmissionProcess>(nodeTransmission, mstp3, inf3, ps3);
    auto tp4 = std::make_shared<TransmissionProcess>(nodeTransmission, mstp4, inf4, ps4);


    geoCOIProb->saveState(1);
    geoCOIProb->setValue(.5);
    EXPECT_TRUE(geoCOI->isDirty());
    EXPECT_TRUE(mstp1->isDirty());
    EXPECT_TRUE(mstp2->isDirty());
    EXPECT_TRUE(mstp3->isDirty());
    EXPECT_TRUE(mstp4->isDirty());
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());

    geoCOIProb->acceptState();
    EXPECT_FALSE(geoCOI->isDirty());
    EXPECT_FALSE(mstp1->isDirty());
    EXPECT_FALSE(mstp2->isDirty());
    EXPECT_FALSE(mstp3->isDirty());
    EXPECT_FALSE(mstp4->isDirty());
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());

    ztmbAssoc->saveState(1);
    ztmbAssoc->setValue(.5);
    EXPECT_TRUE(ztmbCTP->isDirty());
    EXPECT_TRUE(nodeTransmission->isDirty());
    EXPECT_FALSE(mstp1->isDirty());
    EXPECT_FALSE(mstp2->isDirty());
    EXPECT_FALSE(mstp3->isDirty());
    EXPECT_FALSE(mstp4->isDirty());
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    ztmbAssoc->acceptState();
    EXPECT_FALSE(ztmbCTP->isDirty());
    EXPECT_FALSE(mstp1->isDirty());
    EXPECT_FALSE(mstp2->isDirty());
    EXPECT_FALSE(mstp3->isDirty());
    EXPECT_FALSE(mstp4->isDirty());
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());


    ztmbProb->saveState(1);
    ztmbProb->setValue(.5);
    EXPECT_TRUE(ztmbCTP->isDirty());
    EXPECT_TRUE(nodeTransmission->isDirty());
    EXPECT_FALSE(mstp1->isDirty());
    EXPECT_FALSE(mstp2->isDirty());
    EXPECT_FALSE(mstp3->isDirty());
    EXPECT_FALSE(mstp4->isDirty());
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    ztmbProb->acceptState();
    EXPECT_FALSE(ztmbCTP->isDirty());
    EXPECT_FALSE(mstp1->isDirty());
    EXPECT_FALSE(mstp2->isDirty());
    EXPECT_FALSE(mstp3->isDirty());
    EXPECT_FALSE(mstp4->isDirty());
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());

    geoGenProb->saveState(1);
    geoGenProb->setValue(.25);
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    geoGenProb->acceptState();
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());

    afc->alleleFrequencies(as1)->saveState(1);
    afc->alleleFrequencies(as1)->setValue(Simplex({.01, .01, .999996, .01, .01}));
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    afc->alleleFrequencies(as1)->acceptState();
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());


    afc->alleleFrequencies(as2)->saveState(1);
    afc->alleleFrequencies(as2)->setValue(Simplex({.01, .01, .9999999996, .01, .01, .01}));
    EXPECT_TRUE(tp1->isDirty());
    EXPECT_TRUE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    afc->alleleFrequencies(as2)->acceptState();
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());


    infectionOrder->saveState(1);
    infectionOrder->swap(2, 3);
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_TRUE(tp3->isDirty());
    EXPECT_TRUE(tp4->isDirty());
    infectionOrder->acceptState();
    EXPECT_FALSE(tp1->isDirty());
    EXPECT_FALSE(tp2->isDirty());
    EXPECT_FALSE(tp3->isDirty());
    EXPECT_FALSE(tp4->isDirty());
}