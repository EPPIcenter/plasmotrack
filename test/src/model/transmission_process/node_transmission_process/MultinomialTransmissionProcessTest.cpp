//
// Checks that the memoized MultinomialTransmissionProcess returns exactly the same log-likelihoods
// as the implementation it replaced (frozen in MultinomialTransmissionProcessReference.h).
//

#include "gtest/gtest.h"

#include "MultinomialTransmissionProcessReference.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"
#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/parameters/Parameter.h"

#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

#include <boost/random.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

namespace {
    // Same configuration as the model (src/impl/model/Model/config.h).
    constexpr int MAX_ALLELES         = 128;
    constexpr int MAX_COI             = 20;
    constexpr int MAX_PARENTS         = 2;
    constexpr int MAX_STRAINS         = 12;
    constexpr int MAX_PARENT_SET_SIZE = MAX_PARENTS + 1;

    using GeneticsImpl             = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent           = Infection<GeneticsImpl>;
    using AFContainer              = AlleleFrequencyContainer<Simplex, Locus>;
    using COIProbabilityImpl       = ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl   = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AFContainer, InfectionEvent::GenotypeParameterMap, MAX_COI>;
    using ParentSetSizeImpl        = ZTGeometric<MAX_PARENT_SET_SIZE>;
    using Optimized                = MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionImpl, ParentSetSizeImpl, GeneticsImpl>;
    using Reference                = transmission_nets::test::reference::MultinomialTransmissionProcessReference<MAX_PARENTS, MAX_STRAINS, SourceTransmissionImpl, ParentSetSizeImpl>;
    using Engine                   = boost::random::mt19937;

    // allowEmpty: all-zero genotypes exercise the impossible-set paths. Latent parents never carry one in the
    // model, and MultinomialSourceTransmissionProcess exits on the NaN it would produce.
    std::string randomGenotype(Engine& rng, const unsigned int numAlleles, const bool allowEmpty = true) {
        boost::random::uniform_01<> unif;
        std::string bits(numAlleles, '0');
        const double roll = unif(rng);
        if (allowEmpty && roll < 0.03) {
            return bits;
        }
        // Mostly monoclonal, sometimes polyclonal, occasionally more alleles than MAX_STRAINS.
        const unsigned int target = roll < 0.75 ? 1 : (roll < 0.97 ? 2 + static_cast<unsigned int>(unif(rng) * 3) : numAlleles);
        boost::random::uniform_int_distribution<> pick(0, static_cast<int>(numAlleles) - 1);
        for (unsigned int k = 0; k < std::min(target, numAlleles); ++k) {
            bits[pick(rng)] = '1';
        }
        return bits;
    }

    void setLatent(const std::shared_ptr<InfectionEvent>& inf, const std::shared_ptr<Locus>& locus, const std::string& bits) {
        auto param = inf->latentGenotype(locus);
        param->saveState(1);
        param->setValue(GeneticsImpl(bits));
        param->acceptState();
    }

    std::string bitsOf(const std::shared_ptr<InfectionEvent>& inf, const std::shared_ptr<Locus>& locus) {
        return inf->latentGenotype(locus)->value().allelesStr();
    }
}// namespace

TEST(MultinomialTransmissionProcessTest, MatchesReferenceImplementationExactly) {
    Engine rng(20261005);
    boost::random::uniform_01<> unif;

    // Includes a monomorphic locus and one with more alleles than MAX_STRAINS.
    const std::vector<unsigned int> alleleCounts{1, 2, 3, 4, 5, 8, 13, 16};
    std::vector<std::shared_ptr<Locus>> loci;
    auto alleleFreqs = std::make_shared<AFContainer>();
    for (size_t l = 0; l < alleleCounts.size(); ++l) {
        loci.push_back(std::make_shared<Locus>("L" + std::to_string(l), alleleCounts[l]));
        alleleFreqs->addLocus(loci.back());
        alleleFreqs->alleleFrequencies(loci.back())->initializeValue(Simplex(std::vector<double>(alleleCounts[l], 1.0 / alleleCounts[l])));
    }

    // A few founder genotypes, each copied with occasional changes, mimics a set of clonal clusters.
    constexpr int numFounders   = 4;
    constexpr int numInfections = 24;
    std::vector<std::vector<std::string>> founders(numFounders);
    for (auto& founder : founders) {
        for (const auto count : alleleCounts) {
            founder.push_back(randomGenotype(rng, count));
        }
    }

    std::vector<std::shared_ptr<InfectionEvent>> infections;
    std::vector<std::shared_ptr<InfectionEvent>> latentParents;
    for (int i = 0; i < numInfections; ++i) {
        const auto& founder = founders[i % numFounders];
        auto inf            = std::make_shared<InfectionEvent>("inf" + std::to_string(i), 100.0 + i, true);
        for (size_t l = 0; l < loci.size(); ++l) {
            const auto bits = unif(rng) < 0.05 ? randomGenotype(rng, alleleCounts[l]) : founder[l];
            inf->addGenetics(loci[l], bits, bits);
        }
        infections.push_back(inf);

        auto latent = std::make_shared<InfectionEvent>(*inf);
        for (size_t l = 0; l < loci.size(); ++l) {
            if (unif(rng) < 0.1 || latent->latentGenotype(loci[l])->value().totalPositiveCount() == 0) {
                setLatent(latent, loci[l], randomGenotype(rng, alleleCounts[l], false));
            }
        }
        latentParents.push_back(latent);
    }

    auto meanCOI        = std::make_shared<Parameter<double>>(1.5);
    auto coiProb        = std::make_shared<COIProbabilityImpl>(meanCOI);
    auto psProb         = std::make_shared<Parameter<double>>(.9);
    auto psp            = std::make_shared<ParentSetSizeImpl>(psProb);
    auto meanTransmitted = std::make_shared<Parameter<double>>(2.0);

    Optimized optimized(meanTransmitted);
    Reference reference(meanTransmitted);

    std::vector<std::shared_ptr<SourceTransmissionImpl>> sources;
    for (const auto& latent : latentParents) {
        sources.push_back(std::make_shared<SourceTransmissionImpl>(coiProb, alleleFreqs, latent->loci(), latent->latentGenotype()));
    }

    int comparisons = 0;
    int finite      = 0;
    auto expectSame = [&](const double expected, const double actual, const std::string& what) {
        ++comparisons;
        if (std::isfinite(expected)) {
            ++finite;
        }
        // Bit-identical, including -inf.
        EXPECT_EQ(expected, actual) << what;
    };

    auto compareAll = [&](const std::string& round) {
        for (int c = 0; c < numInfections; ++c) {
            const auto& child  = infections[c];
            const auto& latent = latentParents[c];
            const auto& stp    = sources[c];

            expectSame(reference.calculateLogLikelihood(child, latent, stp, psp),
                       optimized.calculateLogLikelihood(child, latent, stp, psp),
                       round + " latent-only child=" + child->id());

            for (int a = 0; a < numInfections; ++a) {
                if (a == c) continue;
                ParentSet<InfectionEvent> single{infections[a]};
                expectSame(reference.calculateLogLikelihood(child, single, psp),
                           optimized.calculateLogLikelihood(child, single, psp),
                           round + " observed child=" + child->id() + " parents=" + infections[a]->id());
                expectSame(reference.calculateLogLikelihood(child, latent, single, stp, psp),
                           optimized.calculateLogLikelihood(child, latent, single, stp, psp),
                           round + " observed+latent child=" + child->id() + " parents=" + infections[a]->id());

                for (int b = a + 1; b < numInfections; b += 3) {
                    if (b == c) continue;
                    ParentSet<InfectionEvent> pair{infections[a], infections[b]};
                    expectSame(reference.calculateLogLikelihood(child, pair, psp),
                               optimized.calculateLogLikelihood(child, pair, psp),
                               round + " observed child=" + child->id() + " parents=" + infections[a]->id() + "," + infections[b]->id());
                    expectSame(reference.calculateLogLikelihood(child, latent, pair, stp, psp),
                               optimized.calculateLogLikelihood(child, latent, pair, stp, psp),
                               round + " observed+latent child=" + child->id() + " parents=" + infections[a]->id() + "," + infections[b]->id());
                }
            }
        }
    };

    compareAll("cold");
    // Second pass is served from the memo.
    compareAll("warm");

    // Change genotypes, so stale memo entries would show up as mismatches, and change a
    // parameter that feeds the strain-count prior rather than the locus terms.
    for (int i = 0; i < numInfections; ++i) {
        const auto& locus = loci[i % loci.size()];
        setLatent(infections[i], locus, randomGenotype(rng, locus->totalAlleles()));
        if (i % 2 == 0) {
            const auto donor = bitsOf(infections[(i + 1) % numInfections], locus);
            setLatent(latentParents[i], locus, donor.find('1') != std::string::npos ? donor : randomGenotype(rng, locus->totalAlleles(), false));
        }
    }
    meanTransmitted->saveState(1);
    meanTransmitted->setValue(3.5);
    meanTransmitted->acceptState();
    compareAll("mutated");

    // Evaluate under a proposed strain-count parameter, then reject it: results cached under the proposal
    // must not be returned once the earlier table is restored.
    meanTransmitted->saveState(1);
    meanTransmitted->setValue(6.0);
    compareAll("proposed");
    meanTransmitted->restoreState(1);
    compareAll("restored");

    EXPECT_GT(comparisons, 1000);
    // Guard against a test that only ever compares -inf.
    EXPECT_GT(finite, comparisons / 10);
}

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
    auto latent = std::make_shared<InfectionEvent>(*child);

    auto alleleFreqs = std::make_shared<AFContainer>();
    alleleFreqs->addLocus(locus);
    alleleFreqs->alleleFrequencies(locus)->initializeValue(Simplex(std::vector<double>{.5, .5}));
    auto coiProb = std::make_shared<COIProbabilityImpl>(std::make_shared<Parameter<double>>(1.5));
    auto stp     = std::make_shared<SourceTransmissionImpl>(coiProb, alleleFreqs, latent->loci(), latent->latentGenotype());
    auto psp     = std::make_shared<ParentSetSizeImpl>(std::make_shared<Parameter<double>>(.9));
    Optimized ntp(std::make_shared<Parameter<double>>(2.0));

    const double twoParentPrior = std::log(psp->value()(2));
    EXPECT_NEAR(ntp.calculateLogLikelihood(child, ParentSet<InfectionEvent>{a, b}, psp), twoParentPrior, 1e-12);
    EXPECT_NEAR(ntp.calculateLogLikelihood(child, latent, ParentSet<InfectionEvent>{a}, stp, psp), stp->value() + twoParentPrior, 1e-12);
    // A single parent has no impossible strain counts, so it is unaffected.
    EXPECT_NEAR(ntp.calculateLogLikelihood(child, ParentSet<InfectionEvent>{a}, psp), std::log(psp->value()(1)), 1e-12);
}
