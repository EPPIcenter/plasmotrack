//
// Created by Maxwell Murphy on 3/10/20.
//
//
#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/distributions/ZTGeometric.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <string>
#include <vector>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::model::transmission_process;

constexpr int MAX_COI     = 12;
constexpr int MAX_ALLELES = 32;

TEST(MultinomialSourceTransmissionProcessTest, BasicTest) {
    using GeneticsImpl             = AllelesBitSet<MAX_ALLELES>;
    using COIProbabilityImpl       = ZTGeometric<MAX_COI>;
    using AlleleFrequencyImpl      = Simplex;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequencyImpl, Locus>;
    using Infection                = Infection<GeneticsImpl, Locus>;

    auto as1 = std::make_shared<Locus>("AS1", 3);
    auto as2 = std::make_shared<Locus>("AS2", 4);

    auto inf1 = std::make_shared<Infection>("inf1", 10.0, false);
    inf1->addGenetics(as1, "001", "001");
    inf1->addGenetics(as2, "0011", "0011");

    auto alleleFreqs = std::make_shared<AlleleFrequencyContainer>();
    alleleFreqs->addLocus(as1);
    alleleFreqs->alleleFrequencies(as1)->initializeValue({.2, .3, .5});
    alleleFreqs->addLocus(as2);
    alleleFreqs->alleleFrequencies(as2)->initializeValue({.1, .2, .3, .4});

    auto coiProb = std::make_shared<Parameter<double>>(.3);
    auto coip    = std::make_shared<COIProbabilityImpl>(coiProb);

    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, Infection::GenotypeParameterMap, MAX_COI> mstp(coip, alleleFreqs, inf1->loci(), inf1->latentGenotype());

    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb->saveState(1);
    EXPECT_FALSE(mstp.isDirty());
    coiProb->setValue(.001);
    EXPECT_TRUE(mstp.isDirty());
    EXPECT_TRUE(coip->isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
    coiProb->restoreState(1);
    EXPECT_FALSE(mstp.isDirty());
    EXPECT_FALSE(coip->isDirty());
    std::cout << "LogLikelihood: " << mstp.value() << std::endl;
}
// The process keeps every locus' row summed and swaps changed rows in, saving and restoring the sums with the rows.
// After many accepted and rejected changes, nested or not, its value must match one computed from scratch.
TEST(MultinomialSourceTransmissionProcessTest, RunningSumsMatchRecomputation) {
    using GeneticsImpl             = AllelesBitSet<MAX_ALLELES>;
    using COIProbabilityImpl       = ZTGeometric<MAX_COI>;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<Simplex, Locus>;
    using Infection                = Infection<GeneticsImpl, Locus>;

    boost::random::mt19937 rng(20261008);
    boost::random::uniform_01<> unif;

    const std::vector<int> alleleCounts{2, 3, 4, 6};
    auto infection   = std::make_shared<Infection>("inf", 10.0, false);
    auto alleleFreqs = std::make_shared<AlleleFrequencyContainer>();
    std::vector<std::shared_ptr<Locus>> loci;
    auto randomGenotype = [&](const int numAlleles) {
        std::string bits(numAlleles, '0');
        boost::random::uniform_int_distribution<> pick(0, numAlleles - 1);
        bits[pick(rng)] = '1';
        if (unif(rng) < 0.5) {
            bits[pick(rng)] = '1';
        }
        return bits;
    };
    auto randomFrequencies = [&](const int numAlleles) {
        std::vector<double> frequencies(numAlleles);
        for (auto& frequency : frequencies) {
            frequency = 0.05 + unif(rng);
        }
        return Simplex(frequencies);
    };
    for (size_t l = 0; l < alleleCounts.size(); ++l) {
        auto locus = std::make_shared<Locus>("L" + std::to_string(l), alleleCounts[l]);
        loci.push_back(locus);
        const auto genotype = randomGenotype(alleleCounts[l]);
        infection->addGenetics(locus, genotype, genotype);
        alleleFreqs->addLocus(locus);
        alleleFreqs->alleleFrequencies(locus)->initializeValue(randomFrequencies(alleleCounts[l]));
    }
    auto coiParameter = std::make_shared<Parameter<double>>(.3);
    auto coiProb      = std::make_shared<COIProbabilityImpl>(coiParameter);

    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, Infection::GenotypeParameterMap, MAX_COI> stp(coiProb, alleleFreqs, infection->loci(), infection->latentGenotype());

    boost::random::uniform_int_distribution<> pickLocus(0, static_cast<int>(loci.size()) - 1);
    // Saves a parameter under `id` and changes it: an allele frequency, a genotype or the COI distribution. Returns
    // functions that restore and accept the change.
    struct Change {
        std::function<void()> restore;
        std::function<void()> accept;
    };
    auto change = [&](const int id) -> Change {
        const auto& locus = loci[pickLocus(rng)];
        auto saved        = [id](const auto& parameter) {
            parameter->saveState(id);
            return Change{.restore = [=]() { parameter->restoreState(id); }, .accept = [=]() { parameter->acceptState(); }};
        };
        const double kind = unif(rng);
        if (kind < 0.45) {
            auto parameter = alleleFreqs->alleleFrequencies(locus);
            auto out       = saved(parameter);
            parameter->setValue(randomFrequencies(static_cast<int>(locus->totalAlleles())));
            return out;
        }
        if (kind < 0.9) {
            auto parameter = infection->latentGenotype(locus);
            auto out       = saved(parameter);
            parameter->setValue(GeneticsImpl(randomGenotype(static_cast<int>(locus->totalAlleles()))));
            return out;
        }
        auto out = saved(coiParameter);
        coiParameter->setValue(0.05 + 0.9 * unif(rng));
        return out;
    };
    auto expectMatches = [&](const std::string& when) {
        const double expected = stp.validate();
        EXPECT_NEAR(stp.value(), expected, 1e-9 * std::max(1.0, std::abs(expected))) << when;
        // Recombined from the running sums, which a missed restore would leave stale.
        stp.setDirty();
        EXPECT_NEAR(stp.value(), expected, 1e-9 * std::max(1.0, std::abs(expected))) << when << ", recombined";
    };

    // Enough swaps that the sums are recomputed from the rows several times.
    for (int k = 0; k < 4000; ++k) {
        const auto when = "step " + std::to_string(k);
        // Sometimes two changes under one saved state, so that one level of the process's undo log holds two rows.
        std::vector<Change> outer{change(1)};
        stp.value();
        if (unif(rng) < 0.3) {
            outer.push_back(change(1));
            stp.value();
        }
        std::optional<Change> inner;
        if (unif(rng) < 0.3) {
            // A nested change, rejected at once or decided with the outer one.
            inner = change(2);
            stp.value();
            if (unif(rng) < 0.5) {
                inner->restore();
                inner.reset();
                expectMatches(when + ", inner restored");
            }
        }
        if (unif(rng) < 0.5) {
            for (const auto& c : outer) {
                c.accept();
            }
            if (inner) {
                inner->accept();
            }
        } else {
            if (inner) {
                inner->restore();
            }
            for (const auto& c : outer) {
                c.restore();
            }
        }
        expectMatches(when);
    }
}
