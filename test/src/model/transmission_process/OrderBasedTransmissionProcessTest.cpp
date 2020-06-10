//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Ordering.h"

#include "core/containers/Infection.h"
#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Locus.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionNoMutation.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

constexpr int MAX_PARENTS = 1;
constexpr int MAX_ALLELES = 32;
constexpr int MAX_COI = 10;
constexpr int MAX_TRANSMISSIONS = 5;

TEST(OrderBasedTransmissionProcessTest, CoreTest) {
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = Infection<GeneticsImpl>;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<Simplex>;
    using Ordering = Ordering<InfectionEvent>;

    using COITransitionProbImpl = ZTMultiplicativeBinomial<MAX_COI>;
    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>;

    using COIProbabilityImpl = ZTGeometric<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEvent, MAX_COI>;
    using TransmissionProcess = OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;

    Locus as1("AS1", 5);
    Locus as2("AS2", 6);

    std::vector<InfectionEvent::LocusGeneticsAssignment> dlas{
            {&as1, GeneticsImpl("11010")},
            {&as2, GeneticsImpl("000011")}
    };

    std::vector<InfectionEvent::LocusGeneticsAssignment> plas{
            {&as1, GeneticsImpl("11010")},
            {&as2, GeneticsImpl("000011")}
    };

    InfectionEvent inf1("inf1", dlas, plas);
    InfectionEvent inf2("inf2", dlas, plas);
    InfectionEvent inf3("inf3", dlas, plas);
    InfectionEvent inf4("inf4", dlas, plas);

    Ordering infectionOrder({&inf1, &inf2, &inf3, &inf4});

    OrderDerivedParentSet ps1(infectionOrder, inf1);
    OrderDerivedParentSet ps2(infectionOrder, inf2);
    OrderDerivedParentSet ps3(infectionOrder, inf3);
    OrderDerivedParentSet ps4(infectionOrder, inf4);

//    std::vector<AlleleFrequencyContainer::LocusAlleleFrequencyAssignment> lfas {
//            {&as1, Simplex(as1.totalAlleles())},
//            {&as2, Simplex(as2.totalAlleles())}
//    };

    AlleleFrequencyContainer afc;
    afc.addLocus(as1);
    afc.addLocus(as2);

    Parameter<double> ztmbProb(.3);
    Parameter<double> ztmbAssoc(1.0);
    COITransitionProbImpl ztmbCTP(ztmbProb, ztmbAssoc);

    Parameter<double> geoGenProb(.8);
    InterTransmissionProbImpl geoGen(geoGenProb);

    NodeTransmissionImpl nodeTransmission(ztmbCTP, geoGen);

    Parameter<double> geoCOIProb(.9);
    COIProbabilityImpl geoCOI(geoCOIProb);

    std::cout << "AFC Addr: " << &afc << std::endl;
    SourceTransmissionImpl mstp1(geoCOI, afc, inf1);
    SourceTransmissionImpl mstp2(geoCOI, afc, inf2);
    SourceTransmissionImpl mstp3(geoCOI, afc, inf3);
    SourceTransmissionImpl mstp4(geoCOI, afc, inf4);

    TransmissionProcess tp1(nodeTransmission, mstp1, inf1, ps1);
    TransmissionProcess tp2(nodeTransmission, mstp2, inf2, ps2);
    TransmissionProcess tp3(nodeTransmission, mstp3, inf3, ps3);
    TransmissionProcess tp4(nodeTransmission, mstp4, inf4, ps4);

    EXPECT_TRUE(tp1.isDirty());
    tp1.value();
    EXPECT_FALSE(tp1.isDirty());

    EXPECT_TRUE(tp2.isDirty());
    tp2.value();
    EXPECT_FALSE(tp2.isDirty());


    geoCOIProb.saveState();
    geoCOIProb.setValue(.5);
    EXPECT_TRUE(geoCOI.isDirty());
    EXPECT_TRUE(mstp1.isDirty());
    EXPECT_TRUE(mstp2.isDirty());
    EXPECT_TRUE(mstp3.isDirty());
    EXPECT_TRUE(mstp4.isDirty());
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());

    geoCOIProb.acceptState();
    EXPECT_FALSE(geoCOI.isDirty());
    EXPECT_FALSE(mstp1.isDirty());
    EXPECT_FALSE(mstp2.isDirty());
    EXPECT_FALSE(mstp3.isDirty());
    EXPECT_FALSE(mstp4.isDirty());
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());

    ztmbAssoc.saveState();
    ztmbAssoc.setValue(.5);
    EXPECT_TRUE(ztmbCTP.isDirty());
    EXPECT_TRUE(nodeTransmission.isDirty());
    EXPECT_FALSE(mstp1.isDirty());
    EXPECT_FALSE(mstp2.isDirty());
    EXPECT_FALSE(mstp3.isDirty());
    EXPECT_FALSE(mstp4.isDirty());
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    ztmbAssoc.acceptState();
    EXPECT_FALSE(ztmbCTP.isDirty());
    EXPECT_FALSE(mstp1.isDirty());
    EXPECT_FALSE(mstp2.isDirty());
    EXPECT_FALSE(mstp3.isDirty());
    EXPECT_FALSE(mstp4.isDirty());
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());


    ztmbProb.saveState();
    ztmbProb.setValue(.5);
    EXPECT_TRUE(ztmbCTP.isDirty());
    EXPECT_TRUE(nodeTransmission.isDirty());
    EXPECT_FALSE(mstp1.isDirty());
    EXPECT_FALSE(mstp2.isDirty());
    EXPECT_FALSE(mstp3.isDirty());
    EXPECT_FALSE(mstp4.isDirty());
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    ztmbProb.acceptState();
    EXPECT_FALSE(ztmbCTP.isDirty());
    EXPECT_FALSE(mstp1.isDirty());
    EXPECT_FALSE(mstp2.isDirty());
    EXPECT_FALSE(mstp3.isDirty());
    EXPECT_FALSE(mstp4.isDirty());
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());

    geoGenProb.saveState();
    geoGenProb.setValue(.25);
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    geoGenProb.acceptState();
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());

    afc.alleleFrequencies(as1).saveState();
    afc.alleleFrequencies(as1).setValue(Simplex({.01, .01, .999996, .01, .01}));
    std::cout << "Allele Freqs address: " << &(afc.alleleFrequencies(as1)) << std::endl;
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    afc.alleleFrequencies(as1).acceptState();
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());


    afc.alleleFrequencies(as2).saveState();
    afc.alleleFrequencies(as2).setValue(Simplex({.01, .01, .9999999996, .01, .01, .01}));
    EXPECT_TRUE(tp1.isDirty());
    EXPECT_TRUE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    afc.alleleFrequencies(as2).acceptState();
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());


    infectionOrder.saveState();
    infectionOrder.swap(2, 3);
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_TRUE(tp3.isDirty());
    EXPECT_TRUE(tp4.isDirty());
    infectionOrder.acceptState();
    EXPECT_FALSE(tp1.isDirty());
    EXPECT_FALSE(tp2.isDirty());
    EXPECT_FALSE(tp3.isDirty());
    EXPECT_FALSE(tp4.isDirty());

}