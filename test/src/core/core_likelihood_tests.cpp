//
// Created by Maxwell Murphy on 12/9/19.
//

#include <optional>
#include <iostream>
#include <boost/container/flat_set.hpp>

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"
#include "core/parameters/Ordering.h"

#include "core/containers/Infection.h"

#include "core/datatypes/Data.h"
#include "core/datatypes/Alleles.h"

#include "core/computation/Accumulator.h"

#include "core/utils/CombinationIndicesGenerator.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/TransmissionProcessLikelihood.h"
#include "model/transmission_process/node_transmission_process/ZTMultiplicativeBinomial.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfection.h"
#include "model/transmission_process/OrderDerivedParentSet.h"
#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/GeometricGenerationProbability.h"


using GeneticsImpl = AllelesBitSet<16>;
using DoubleParameter = Parameter<double>;

TEST(CoreLikelihoodTest, LikelihoodTest) {

    using AlleleCounterAccumulator = Accumulator<AlleleCounter<GeneticsImpl>, AlleleCounts>;

    Data<GeneticsImpl> d1("011010");
    Data<GeneticsImpl> d2("101101");
    Parameter<GeneticsImpl> k1("011010");
    Parameter<GeneticsImpl> k2("101101");

    std::cout << "Created Parameters" << std::endl;
    AlleleCounter ac1(k1, d1);
    AlleleCounter ac2(k2, d2);

    std::cout << "Created Allele Counters" << std::endl;

    AlleleCounterAccumulator acc;
    acc.addTarget(ac1);
    acc.addTarget(ac2);

    std::cout << "Added Allele Counter Targets" << std::endl;

    std::cout << ac1.value() << std::endl;
    std::cout << ac2.value() << std::endl;
    std::cout << "Accumulator:" << std::endl;
    std::cout << acc.value() << std::endl;

    DoubleParameter fpr(.05);
    DoubleParameter fnr(.05);

    ObservationProcessLikelihood op(acc, fpr, fnr);

    Accumulator<PartialLikelihood, float> llik;
    llik.addTarget(op);

    std::cout << llik.value() << std::endl;

    fpr.saveState();
    std::cout << "setting fpr value" << std::endl;
    fpr.setValue(.1);
    std::cout << "set fpr value" << std::endl;
    std::cout << llik.value() << std::endl;

    fpr.restoreState();
    std::cout << llik.value() << std::endl;

    fnr.saveState();
    fnr.setValue(.99);
    std::cout << llik.value() << std::endl;
    fnr.acceptState();

    int a = 3;
    int b = 4;
    int c = 6;
    int d = 7;

    Ordering<int> op2;
    op2.addElements({&a, &b, &c, &d});
    OrderDerivedParentSet ps(op2, c);

    std::cout << "Printing Parent Set" << std::endl;
    ps.printSet();
    std::cout << op2 << std::endl;
    op2.swap(0, 3);
    std::cout << "Printing Parent Set" << std::endl;
    ps.printSet();
    std::cout << op2 << std::endl;
    op2.swap(2, 3);
    std::cout << "Printing Parent Set" << std::endl;
    ps.printSet();
    std::cout << op2 << std::endl;
    op2.swap(3, 0);
    std::cout << "Printing Parent Set" << std::endl;
    ps.printSet();
    std::cout << op2 << std::endl;

    DoubleParameter tcp(.85);
    DoubleParameter tca(.95);
    ZTMultiplicativeBinomial<10> ztmb(tcp, tca);
    std::cout << ztmb.value()(5, 3) << std::endl;

    tcp.saveState();
    tca.saveState();
    tcp.setValue(.95);
    tca.setValue(2.0);
    std::cout << ztmb.value()(5, 3) << std::endl;
    tcp.restoreState();
    tca.restoreState();
    std::cout << ztmb.value() << std::endl;

    DoubleParameter gp_prob(.6);
    GeometricGenerationProbability<25> gp(gp_prob);

    NoSuperInfection<10, 25, ZTMultiplicativeBinomial, GeometricGenerationProbability> tp(ztmb, gp);

    std::cout << "Log Probability: " << std::endl;
    std::cout << tp.value() << std::endl;

    Locus as1("AS1", 6);

    std::vector<LocusGeneticsAssignment<GeneticsImpl, Locus>> dlas{{&as1, GeneticsImpl("011010")}};
    std::vector<LocusGeneticsAssignment<GeneticsImpl, Locus>> plas{{&as1, GeneticsImpl("011010")}};
    Infection<GeneticsImpl> inf1(dlas, plas);
    Infection<GeneticsImpl> inf2(dlas, plas);
    Infection<GeneticsImpl> inf3(dlas, plas);
    Infection<GeneticsImpl> inf4(dlas, plas);

    ParentSet<Infection<GeneticsImpl>> ps1{&inf1};

    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    tcp.saveState();
    tcp.setValue(.05);
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    tcp.restoreState();
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;

    inf1.latentGenotype(&as1).saveState();
    std::cout << "Saved State" << std::endl;
    inf1.latentGenotype(&as1).setValue(GeneticsImpl("111111"));
    std::cout << "Set Value" << std::endl;
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    inf1.latentGenotype(&as1).restoreState();
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;


    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;
    gp_prob.saveState();
    gp_prob.setValue(.1);
    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;
    gp_prob.restoreState();
    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;


}
