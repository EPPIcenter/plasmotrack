//
// Created by Maxwell Murphy on 12/9/19.
//

#include <optional>
#include <iostream>
#include <boost/container/flat_set.hpp>

#include "gtest/gtest.h"

#include "core/parameter/Parameter.h"
#include "core/parameter/Ordering.h"

#include "core/datatypes/Data.h"
#include "core/datatypes/Alleles.h"

#include "core/computation/Accumulator.h"

#include "core/computation/observation_process/AlleleCounter.h"
#include "core/computation/observation_process/ObservationProcessLikelihood.h"
#include "core/computation/transmission_process/TransmissionProcessLikelihood.h"
#include "core/computation/transmission_process/ZTMultiplicativeBinomial.h"
#include "core/computation/transmission_process/TransmissionProcess.h"
#include "core/computation/transmission_process/OrderDerivedParentSet.h"
#include "core/computation/transmission_process/OrderBasedTransmissionProcess.h"


using GeneticsImpl = AllelesBitSet<16>;

TEST(CoreLikelihoodTest, LikelihoodTest) {

    using AlleleCounterAccumulator = Accumulator<AlleleCounter<GeneticsImpl>, AlleleCounts>;

    Data d1("d1", GeneticsImpl("011010"));
    Data d2("d2", GeneticsImpl("101101"));
    Parameter k1("a1", d1.value());
    Parameter k2("a1", d2.value());

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

    Parameter fpr("fpr", .05);
    Parameter fnr("fnr", .05);

    ObservationProcessLikelihood op("ob", acc, fpr, fnr);

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

//    auto a = new int{3};
//    auto b = new int{4};
//    auto c = new int{5};
//    auto d = new int{6};

    int a = 3;
    int b = 4;
    int c = 6;
    int d = 7;

    Ordering<int> op2("op2");
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

    Parameter tcp("tcp", .85);
    Parameter tca("tca", .95);
    ZTMultiplicativeBinomial<20> ztmb(tcp, tca);
    std::cout << ztmb.value()(5, 3) << std::endl;

    tcp.saveState();
    tca.saveState();
    tcp.setValue(.95);
    tca.setValue(2.0);
    std::cout << ztmb.value()(5, 3) << std::endl;
    tcp.restoreState();
    tca.restoreState();
    std::cout << ztmb.value() << std::endl;

    TransmissionProcess<3, 20> tp(ztmb);

    std::cout << "Probability: " << std::endl;
    std::cout << tp.value() << std::endl;

}
