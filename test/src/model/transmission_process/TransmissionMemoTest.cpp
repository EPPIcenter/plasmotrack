//
// The Transmission memo keeps kernel results by genome and frequency identity: Infections with the same genotypes
// share ids and sums, a Locus's frequencies get an id that follows their values, the strain-count table follows its
// parameter, and the latent locus tables at the frequencies a term last summed with stay available.
//

#include "gtest/gtest.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Genotype.h"
#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "model/transmission_process/TransmissionMemo.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;
using namespace transmission_nets::model::transmission_process;

namespace {
    constexpr int MAX_PARENTS = 2;
    constexpr int MAX_STRAINS = 12;
    constexpr int MAX_COI     = 20;
    using InfectionEvent = transmission_nets::core::containers::Infection<>;
    using Frequencies = AlleleFrequencyContainer<Simplex, Locus>;
    using Memo        = TransmissionMemo<MAX_PARENTS, MAX_STRAINS, MAX_COI, InfectionEvent, Locus, Frequencies>;
    using Kernel      = Memo::Kernel;
    using Variant     = Memo::Variant;

    struct Fixture {
        std::vector<std::shared_ptr<Locus>> loci{std::make_shared<Locus>("L0", 3), std::make_shared<Locus>("L1", 4)};
        std::shared_ptr<Frequencies> frequencies = std::make_shared<Frequencies>();
        std::shared_ptr<Parameter<double>> meanStrains = std::make_shared<Parameter<double>>(1.7);
        Memo memo{meanStrains, frequencies};

        Fixture() {
            frequencies->addLocus(loci[0]);
            frequencies->addLocus(loci[1]);
            frequencies->alleleFrequencies(loci[0])->initializeValue(Simplex(std::vector<double>{0.5, 0.3, 0.2}));
            frequencies->alleleFrequencies(loci[1])->initializeValue(Simplex(std::vector<double>{0.4, 0.3, 0.2, 0.1}));
        }
        // An Infection carrying `bits0` and `bits1` at the two loci, observed and latent alike.
        std::shared_ptr<InfectionEvent> infection(const std::string& id, const std::string& bits0, const std::string& bits1) const {
            auto out = std::make_shared<InfectionEvent>(id, 100.0, true);
            out->addGenetics(loci[0], bits0, bits0);
            out->addGenetics(loci[1], bits1, bits1);
            return out;
        }
        static void setLatent(const std::shared_ptr<InfectionEvent>& inf, const std::shared_ptr<Locus>& locus, const std::string& bits) {
            auto genotype = inf->latentGenotype(locus);
            genotype->saveState(1);
            genotype->setValue(Genotype(bits));
            genotype->acceptState();
        }
        void setFrequencies(const std::shared_ptr<Locus>& locus, const std::vector<double>& values) const {
            auto parameter = frequencies->alleleFrequencies(locus);
            parameter->saveState(1);
            parameter->setValue(Simplex(values));
            parameter->acceptState();
        }
    };
}// namespace

TEST(TransmissionMemoTest, InfectionsWithTheSameLatentGenotypesShareAGenomeIdUntilOneChanges) {
    Fixture f;
    const auto a = f.infection("a", "100", "1000");
    const auto b = f.infection("b", "100", "1000");
    const auto c = f.infection("c", "010", "1000");
    EXPECT_EQ(f.memo.genomeOf(a), f.memo.genomeOf(b));
    EXPECT_NE(f.memo.genomeOf(a), f.memo.genomeOf(c));

    Fixture::setLatent(b, f.loci[0], "010");
    EXPECT_EQ(f.memo.genomeOf(b), f.memo.genomeOf(c)) << "the id follows the genotypes, not the Infection";
    Fixture::setLatent(b, f.loci[0], "100");
    EXPECT_EQ(f.memo.genomeOf(b), f.memo.genomeOf(a));
}

TEST(TransmissionMemoTest, AFrequencyIdFollowsTheFrequenciesValues) {
    Fixture f;
    const auto first = f.memo.frequencyId(f.loci[0]);
    EXPECT_EQ(f.memo.frequencyId(f.loci[0]), first);
    f.setFrequencies(f.loci[0], {0.2, 0.3, 0.5});
    const auto second = f.memo.frequencyId(f.loci[0]);
    EXPECT_NE(second, first);
    f.setFrequencies(f.loci[0], {0.5, 0.3, 0.2});
    EXPECT_EQ(f.memo.frequencyId(f.loci[0]), first) << "the same values again get the same id";
    EXPECT_EQ(f.memo.generation(), 0u);
}

TEST(TransmissionMemoTest, LocusTermsAndParentSetSumsAreSharedByIdenticalGenotypes) {
    Fixture f;
    f.memo.beginBatch();
    const auto child  = f.infection("child", "110", "1100");
    const auto twin   = f.infection("twin", "110", "1100");
    const auto parent = f.infection("parent", "110", "1100");
    const Genotype childAt0 = child->latentGenotype(f.loci[0])->value(), parentAt0 = parent->latentGenotype(f.loci[0])->value();
    const std::vector<const Genotype*> parents{&parentAt0};
    EXPECT_EQ(&f.memo.locusTerm(Variant::ObservedParents, childAt0, parents, nullptr), &f.memo.locusTerm(Variant::ObservedParents, childAt0, parents, nullptr));

    const ParentSet<InfectionEvent> parentSet{parent};
    Memo::ResolvedGenomes forChild{.childLocusOrder = f.memo.locusOrderOf(child), .child = f.memo.genomeOf(child), .parents = {f.memo.genomeOf(parent), 0}};
    Memo::ResolvedGenomes forTwin{.childLocusOrder = f.memo.locusOrderOf(twin), .child = f.memo.genomeOf(twin), .parents = {f.memo.genomeOf(parent), 0}};
    EXPECT_EQ(forChild.child, forTwin.child);
    EXPECT_EQ(&f.memo.parentSetSums(Variant::ObservedParents, child, &parentSet, nullptr, forChild), &f.memo.parentSetSums(Variant::ObservedParents, twin, &parentSet, nullptr, forTwin)) << "one entry for both children";
}

TEST(TransmissionMemoTest, TheStrainCountTableFollowsItsParameterThroughChangeAndRestore) {
    Fixture f;
    EXPECT_EQ(f.memo.strainCountTable(), Kernel::strainCountTable(1.7));
    f.meanStrains->saveState(1);
    f.meanStrains->setValue(2.5);
    EXPECT_EQ(f.memo.strainCountTable(), Kernel::strainCountTable(2.5));
    f.meanStrains->restoreState(1);
    EXPECT_EQ(f.memo.strainCountTable(), Kernel::strainCountTable(1.7));
}

TEST(TransmissionMemoTest, LatentLocusTablesAtEarlierFrequenciesStayAvailable) {
    Fixture f;
    f.memo.beginBatch();
    const Genotype child("110");
    const auto before = f.memo.frequencyId(f.loci[0]);
    const auto atBefore = f.memo.locusTable(before, child, {}).terms;
    f.setFrequencies(f.loci[0], {0.2, 0.3, 0.5});
    const auto after = f.memo.frequencyId(f.loci[0]);
    const auto atAfter = f.memo.locusTable(after, child, {}).terms;
    EXPECT_NE(atBefore, atAfter);
    // An allele-frequency proposal needs both: the table at the frequencies last summed with, to swap out, and the one
    // at the proposed frequencies, to swap in.
    EXPECT_EQ(f.memo.locusTable(before, child, {}).terms, atBefore);
    EXPECT_EQ(f.memo.locusTable(after, child, {}).terms, atAfter);
}

TEST(TransmissionMemoTest, ThePlainObservedParentsTermIsTheKernelsSumOverLoci) {
    Fixture f;
    const auto child  = f.infection("child", "110", "1100");
    const auto parent = f.infection("parent", "100", "1100");
    Kernel kernel;
    Kernel::StrainSums sums{};
    for (const auto& locus : f.loci) {
        const Genotype parentAt = parent->latentGenotype(locus)->value();
        const std::vector<const Genotype*> parents{&parentAt};
        sums.add(kernel.locusTerm(Variant::ObservedParents, child->latentGenotype(locus)->value(), parents, nullptr));
    }
    const double expected = Kernel::strainSummedLogLikelihood(sums, 1, Kernel::strainCountTable(1.7));
    EXPECT_DOUBLE_EQ(f.memo.observedParentsLogLikelihood(child, ParentSet<InfectionEvent>{parent}), expected);
}
