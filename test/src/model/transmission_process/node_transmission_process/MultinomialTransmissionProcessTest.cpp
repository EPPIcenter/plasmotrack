//
// The strain-count sum of the node transmission term: counts below the number of parents are impossible and must not
// contribute.
//

#include "gtest/gtest.h"

#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Alleles.h"
#include "core/distributions/ZTGeometric.h"
#include "core/parameters/Parameter.h"

#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h"

#include <cmath>
#include <memory>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

namespace {
    // Same configuration as the model (src/impl/model/Model/config.h).
    constexpr int MAX_ALLELES         = 128;
    constexpr int MAX_PARENTS         = 2;
    constexpr int MAX_STRAINS         = 12;
    constexpr int MAX_PARENT_SET_SIZE = MAX_PARENTS + 1;

    using GeneticsImpl      = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent    = Infection<GeneticsImpl>;
    using ParentSetSizeImpl = ZTGeometric<MAX_PARENT_SET_SIZE>;
    using NodeTerm          = MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizeImpl, GeneticsImpl>;
}// namespace

TEST(MultinomialTransmissionProcessTest, ExcludesStrainCountsBelowNumberOfParents) {
    // Every strain carries the child's only allele, so each locus term is log(1) = 0 and the likelihood reduces to
    // the strain-count distribution summed over possible counts (= 1) plus the parent-set-size prior. Counting the
    // impossible single-strain case for two parents used to add exactly log(2).
    auto locus = std::make_shared<Locus>("L0", 2);
    auto child = std::make_shared<InfectionEvent>("child", 100.0, true);
    auto a     = std::make_shared<InfectionEvent>("a", 90.0, true);
    auto b     = std::make_shared<InfectionEvent>("b", 80.0, true);
    for (const auto& inf : {child, a, b}) {
        inf->addGenetics(locus, "10", "10");
    }

    auto psp = std::make_shared<ParentSetSizeImpl>(std::make_shared<Parameter<double>>(.9));
    NodeTerm ntp(std::make_shared<Parameter<double>>(2.0));

    const double twoParentPrior = std::log(psp->value()(2));
    EXPECT_NEAR(ntp.calculateLogLikelihood(child, ParentSet<InfectionEvent>{a, b}, psp), twoParentPrior, 1e-12);
    // A single parent has no impossible strain counts, so it is unaffected.
    EXPECT_NEAR(ntp.calculateLogLikelihood(child, ParentSet<InfectionEvent>{a}, psp), std::log(psp->value()(1)), 1e-12);
}
