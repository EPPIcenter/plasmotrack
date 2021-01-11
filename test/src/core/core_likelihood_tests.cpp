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

#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/general/ContinuousRandomWalk.h"
#include "core/samplers/general/DiscreteRandomWalk.h"
#include "core/samplers/general/SALTSampler.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "core/computation/OrderDerivedParentSet.h"
#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"
#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/TransmissionProcessLikelihood.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionNoMutation.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"


using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::utils;
using namespace transmission_nets::core::samplers;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::observation_process;
using namespace transmission_nets::model::transmission_process;

using GeneticsImpl = AllelesBitSet<16>;
using DoubleParameter = Parameter<double>;

TEST(CoreLikelihoodTest, LikelihoodTest) {

    using AlleleCounterAccumulator = Accumulator<AlleleCounter<GeneticsImpl>, AlleleCounts>;
    using Infection = Infection<GeneticsImpl, Locus>;

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

    fpr.saveState("state1");
    std::cout << "setting fpr value" << std::endl;
    fpr.setValue(.1);
    std::cout << "set fpr value" << std::endl;
    std::cout << llik.value() << std::endl;

    fpr.restoreState("state1");
    std::cout << llik.value() << std::endl;

    fnr.saveState("state1");
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
    ps.serialize();
    std::cout << op2 << std::endl;
    op2.swap(0, 3);
    std::cout << "Printing Parent Set" << std::endl;
    ps.serialize();
    std::cout << op2 << std::endl;
    op2.swap(2, 3);
    std::cout << "Printing Parent Set" << std::endl;
    ps.serialize();
    std::cout << op2 << std::endl;
    op2.swap(3, 0);
    std::cout << "Printing Parent Set" << std::endl;
    ps.serialize();
    std::cout << op2 << std::endl;

    DoubleParameter tcp(.85);
    DoubleParameter tca(.95);
    ZTMultiplicativeBinomial<10> ztmb(tcp, tca);
    std::cout << ztmb.value()(5, 3) << std::endl;

    tcp.saveState("state1");
    tca.saveState("state1");
    tcp.setValue(.95);
    tca.setValue(2.0);
    std::cout << ztmb.value()(5, 3) << std::endl;
    tcp.restoreState("state1");
    tca.restoreState("state1");
    std::cout << ztmb.value() << std::endl;

    DoubleParameter gp_prob(.6);
    ZTGeometric<25> gp(gp_prob);

    NoSuperInfectionNoMutation<10, 25, ZTMultiplicativeBinomial<10>, ZTGeometric<25>> tp(ztmb, gp);

    std::cout << "Log Probability: " << std::endl;
    std::cout << tp.value() << std::endl;

    Locus as1("AS1", 6);

    std::vector<Infection::LocusGeneticsAssignment> dlas{{&as1, GeneticsImpl("011010")}};
    std::vector<Infection::LocusGeneticsAssignment> plas{{&as1, GeneticsImpl("011010")}};
    Infection inf1("inf1", dlas, plas);
    Infection inf2("inf2", dlas, plas);
    Infection inf3("inf3", dlas, plas);
    Infection inf4("inf4", dlas, plas);

    ParentSet<Infection> ps1{&inf1};

    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    tcp.saveState("state1");
    tcp.setValue(.05);
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    tcp.restoreState("state1");
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;

    inf1.latentGenotype(&as1).saveState("state1");
    std::cout << "Saved State" << std::endl;
    inf1.latentGenotype(&as1).setValue(GeneticsImpl("111111"));
    std::cout << "Set Value" << std::endl;
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;
    inf1.latentGenotype(&as1).restoreState("state1");
    std::cout << "Parent Set: " << tp.calculateLikelihood(inf2, ps1) << std::endl;


    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;
    gp_prob.saveState("state1");
    gp_prob.setValue(.1);
    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;
    gp_prob.restoreState("state1");
    std::cout << "Geometric Prob: " << gp.value() << std::endl;
    std::cout << "Geometric Prob Sum: " << gp.value().sum() << std::endl;


}
