//
// The transmission leak as a sampled Run scalar (--leak-prior), and the shared and hierarchical error models.
//
// - Cache consistency: after sweeps of the full Sampler set, the leak's walk included, and forced leak proposals
//   accepted and rejected between genotype and frequency changes, the cached log posterior equals that of a Model
//   built afresh on a copy of the State. A rejected leak step brings the memo's results at the old leak back; this
//   catches a term whose restored cache and the memo's results disagree.
// - Stationarity: the leak's walk alone, everything else fixed, has the conditional posterior of the leak as its
//   stationary distribution, computed exactly on a grid.
//

#include "gtest/gtest.h"

#include "Studies.h"
#include "core/samplers/Proposal.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "impl/model/Model/EnsembleCheckpoint.h"
#include "impl/model/Model/ErrorRatePrior.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/SamplerSet.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using namespace transmission_nets::test::studies;

namespace {
    // Clonal clusters with discordant alleles, so that the leak matters, and missing genotypes.
    Study leakStudy(boost::random::mt19937& rng, const int numNodes = 10) {
        const std::vector<int> alleleCounts{2, 3, 4, 5, 6};
        constexpr int numClusters = 2;
        boost::random::uniform_01<> unif;
        auto randomGenotype = [&](const int numAlleles) {
            std::string bits(numAlleles, '0');
            boost::random::uniform_int_distribution<> pick(0, numAlleles - 1);
            bits[pick(rng)] = '1';
            if (unif(rng) < 0.2) {
                bits[pick(rng)] = '1';
            }
            return bits;
        };
        std::vector<std::vector<std::string>> founders(numClusters);
        for (auto& founder : founders) {
            for (const auto count : alleleCounts) {
                founder.push_back(randomGenotype(count));
            }
        }
        Study study{.symptomaticInfectionDurations = flatDurations(), .asymptomaticInfectionDurations = flatDurations()};
        for (size_t l = 0; l < alleleCounts.size(); ++l) {
            study.loci.push_back(uniformLocus("L" + std::to_string(l), static_cast<unsigned int>(alleleCounts[l])));
        }
        for (int n = 0; n < numNodes; ++n) {
            Study::Infection infection{.id = "n" + std::to_string(n), .observationTime = 400.0 + 15 * n, .symptomatic = n % 3 != 0};
            for (size_t l = 0; l < alleleCounts.size(); ++l) {
                std::string genotype = unif(rng) < 0.2 ? randomGenotype(alleleCounts[l]) : founders[n % numClusters][l];
                if (unif(rng) < 0.05) {
                    genotype = std::string(alleleCounts[l], '0');
                }
                infection.observedGenotypes.push_back({"L" + std::to_string(l), genotype});
            }
            for (int p = 0; p < numNodes; ++p) {
                if (p != n) {
                    infection.permittedParents.push_back("n" + std::to_string(p));
                }
            }
            study.infections.push_back(std::move(infection));
        }
        return study;
    }

    // The log posterior of a Model built afresh on a copy of `state`.
    double freshValue(const Study& study, const ModelChoice& choice, const State& state) {
        auto rng   = std::make_shared<boost::random::mt19937>(1);
        auto fresh = std::make_shared<State>(study, rng, choice);
        checkpoint::restoreState(checkpoint::stateLayout(state), checkpoint::captureState(state), *fresh);
        Model model(fresh);
        return model.value();
    }

    void checkConsistency(const ModelChoice& choice, const unsigned int seed) {
        auto rng         = std::make_shared<boost::random::mt19937>(seed);
        const auto study = leakStudy(*rng);
        auto state       = std::make_shared<State>(study, rng, choice);
        auto model       = std::make_shared<Model>(state);
        SamplerSet<Model> samplers(state, model, rng);

        boost::random::uniform_01<> unif;
        boost::random::normal_distribution<> normal;
        transmission_nets::core::samplers::SamplerCounts counts{};
        const auto perturb = [&](const auto& param, const auto& proposed, const bool accept) {
            auto step = transmission_nets::core::samplers::Proposal(*model, *rng, counts);
            step.change(*param, proposed);
            if (accept and std::isfinite(step.value())) {
                step.accept();
            } else {
                step.reject();
            }
            model->value();
        };
        boost::random::uniform_int_distribution<> pickNode(0, static_cast<int>(state->infections.size()) - 1);
        for (int round = 0; round < 12; ++round) {
            samplers.step();
            for (int k = 0; k < 6; ++k) {
                // A leak step, accepted or rejected, then a genotype change and a frequency change.
                if (state->leakPrior) {
                    const double leak = std::clamp(state->leakSize->value() * std::exp(normal(*rng)), state->leakPrior->lower() * 1.01, state->leakPrior->upper() * 0.99);
                    perturb(state->leakSize, leak, unif(*rng) < 0.5);
                }
                const auto& infection = state->infections[pickNode(*rng)];
                const auto& [label, locus] = *std::next(state->loci.begin(), static_cast<long>(unif(*rng) * static_cast<double>(state->loci.size())));
                auto genotype = infection->latentGenotype(locus)->value();
                genotype.flip(static_cast<size_t>(unif(*rng) * genotype.totalAlleles()));
                if (genotype.totalPositiveCount() > 0 and genotype.totalPositiveCount() < MAX_COI) {
                    perturb(infection->latentGenotype(locus), genotype, unif(*rng) < 0.5);
                }
                auto frequencies = state->alleleFrequencies->alleleFrequencies(locus)->value().frequencies();
                for (auto& f : frequencies) {
                    f *= std::exp(0.2 * normal(*rng));
                }
                double total = 0;
                for (const auto f : frequencies) total += f;
                for (auto& f : frequencies) f /= total;
                perturb(state->alleleFrequencies->alleleFrequencies(locus), AlleleFrequencyImpl::exactly(frequencies), unif(*rng) < 0.5);
            }
            ASSERT_NEAR(model->value(), freshValue(study, choice, *state), 1e-6 * std::max(1.0, std::abs(model->value()))) << "round " << round << ", leak " << state->leakSize->value();
        }
    }
}// namespace

TEST(SampledLeakTest, LogUniformPriorLeakCachesMatchAFreshModel) {
    checkConsistency({.transmissionLeak = 0.01, .leakPrior = "loguniform:1e-6,0.2"}, 20261031);
}

TEST(SampledLeakTest, ExponentialPriorLeakCachesMatchAFreshModel) {
    checkConsistency({.transmissionLeak = 0.01, .leakPrior = "exponential:50"}, 20261032);
}

TEST(SampledLeakTest, SharedAndHierarchicalErrorModelsMatchAFreshModel) {
    checkConsistency({.transmissionLeak = 0.01, .leakPrior = "loguniform:1e-6,0.2", .errorModel = "shared"}, 20261033);
    checkConsistency({.errorModel = "hierarchical"}, 20261034);
}

TEST(SampledLeakTest, TheLeakIsARunScalarOnlyWithAPrior) {
    auto rng         = std::make_shared<boost::random::mt19937>(3);
    const auto study = leakStudy(*rng, 4);
    const auto names = [&](const ModelChoice& choice) { return State(study, rng, choice).scalarNames(); };
    EXPECT_EQ(std::ranges::count(names({.transmissionLeak = 0.002}), "transmission_leak"), 0);
    EXPECT_EQ(std::ranges::count(names({.transmissionLeak = 0.002, .leakPrior = "loguniform"}), "transmission_leak"), 1);
    EXPECT_EQ(std::ranges::count(names({.errorModel = "hierarchical"}), "efn_shape"), 1);
    EXPECT_THROW(parseLeakPrior("loguniform:0.1,0.01"), std::invalid_argument);
    EXPECT_THROW(parseLeakPrior("gamma:2"), std::invalid_argument);
    EXPECT_THROW(parseErrorModel("per-locus"), std::invalid_argument);
    // One walk on each shared rate.
    State shared(study, rng, ModelChoice{.errorModel = "shared"});
    auto sharedState = std::make_shared<State>(shared);
    auto model       = std::make_shared<Model>(sharedState);
    SamplerSet<Model> samplers(sharedState, model, rng);
    EXPECT_EQ(std::ranges::count_if(samplers.ids(), [](const auto& id) { return id.starts_with("ConstrainedContinuousRandomWalk/expected_false_positives/"); }), 1);
}

// The leak's walk alone, with the rest of the State fixed, against its exact conditional on a grid of log ε.
TEST(SampledLeakTest, TheLeakWalkSamplesItsConditionalPosterior) {
    auto rng         = std::make_shared<boost::random::mt19937>(77);
    const auto study = leakStudy(*rng, 8);
    const ModelChoice choice{.transmissionLeak = 0.01, .leakPrior = "loguniform:1e-4,0.5"};
    auto state = std::make_shared<State>(study, rng, choice);
    auto model = std::make_shared<Model>(state);
    {
        // Somewhere typical for the rest of the State.
        SamplerSet<Model> samplers(state, model, rng);
        for (int s = 0; s < 20; ++s) {
            samplers.step();
        }
    }
    const auto& prior = *state->leakPrior;
    transmission_nets::core::samplers::ConstrainedContinuousRandomWalk<Model> walk(state->leakSize, model, prior.lower(), prior.upper(), rng, 1.5);

    // The exact conditional of u = log ε on a grid, from the Model at each ε (its prior term included).
    constexpr int gridSize = 400;
    const double lo = std::log(prior.lower()), hi = std::log(prior.upper());
    std::vector<double> u(gridSize), logDensity(gridSize);
    transmission_nets::core::samplers::SamplerCounts counts{};
    for (int g = 0; g < gridSize; ++g) {
        u[g] = lo + (hi - lo) * (g + 0.5) / gridSize;
        auto step = transmission_nets::core::samplers::Proposal(*model, *rng, counts);
        step.change(*state->leakSize, std::exp(u[g]));
        logDensity[g] = step.value() + u[g];// the Jacobian of ε = exp(u)
        step.accept();
    }
    const double maxLog = *std::ranges::max_element(logDensity);
    std::vector<double> cdf(gridSize);
    double total = 0;
    for (int g = 0; g < gridSize; ++g) {
        total += std::exp(logDensity[g] - maxLog);
        cdf[g] = total;
    }
    for (auto& c : cdf) c /= total;

    constexpr int draws = 25000;
    std::vector<double> sampled;
    sampled.reserve(draws);
    for (int d = 0; d < draws; ++d) {
        walk.update();
        sampled.push_back(std::log(state->leakSize->value()));
    }
    std::ranges::sort(sampled);
    double maxGap = 0;
    for (int g = 0; g < gridSize; g += 10) {
        const double edge      = lo + (hi - lo) * (g + 1.0) / gridSize;
        const double empirical = static_cast<double>(std::ranges::upper_bound(sampled, edge) - sampled.begin()) / draws;
        maxGap                 = std::max(maxGap, std::abs(empirical - cdf[g]));
    }
    EXPECT_LT(maxGap, 0.04) << "largest CDF difference between the walk and the exact conditional";
}

// Controls for the harness: the fixed leak and no leak, as the Model always had them.
TEST(SampledLeakTest, AFixedLeakAndNoLeakMatchAFreshModel) {
    checkConsistency({}, 20261035);
    checkConsistency({.transmissionLeak = 0.01}, 20261036);
}

// Prior recovery: with the likelihood switched off (temperature 0) the whole Sampler set must sample the prior. The
// leak's log-uniform prior on [1e-4, 0.1] is uniform in log10 ε on [-4, -1]; the hierarchical means have Gamma(2, 0.01)
// and Gamma(2, 0.02) priors, means 0.02 and 0.04.
TEST(SampledLeakTest, WithoutTheLikelihoodTheSamplerSetRecoversThePrior) {
    auto rng         = std::make_shared<boost::random::mt19937>(4242);
    const auto study = leakStudy(*rng, 4);
    auto state       = std::make_shared<State>(study, rng, ModelChoice{.transmissionLeak = 0.01, .leakPrior = "loguniform:1e-4,0.1", .errorModel = "hierarchical"});
    auto model       = std::make_shared<Model>(state, 0.0);
    SamplerSet<Model> samplers(state, model, rng);
    std::vector<double> leak, fpMean, fnMean;
    for (int it = 0; it < 3000; ++it) {
        samplers.step(it < 200);
        if (it < 500) {
            continue;
        }
        leak.push_back(std::log10(state->leakSize->value()));
        fpMean.push_back(state->efpMean->value());
        fnMean.push_back(state->efnMean->value());
    }
    // The draws' mean against the prior's, within 4 batch-means standard errors (20 batches).
    const auto expectMean = [](const std::vector<double>& x, const double expected, const char* what) {
        constexpr int batches = 20;
        const size_t size     = x.size() / batches;
        std::vector<double> means(batches);
        for (int b = 0; b < batches; ++b) {
            for (size_t k = 0; k < size; ++k) means[b] += x[b * size + k] / static_cast<double>(size);
        }
        double mean = 0, var = 0;
        for (const double m : means) mean += m / batches;
        for (const double m : means) var += (m - mean) * (m - mean) / (batches - 1);
        EXPECT_LT(std::abs(mean - expected), 4 * std::sqrt(var / batches)) << what << ": mean " << mean << ", prior mean " << expected;
    };
    expectMean(fpMean, 0.02, "efp_mean");
    expectMean(fnMean, 0.04, "efn_mean");
    std::ranges::sort(leak);
    double maxGap = 0;
    for (size_t k = 0; k < leak.size(); ++k) {
        maxGap = std::max(maxGap, std::abs(static_cast<double>(k + 1) / static_cast<double>(leak.size()) - (leak[k] + 4) / 3));
    }
    EXPECT_LT(maxGap, 0.1) << "log10 leak against uniform on [-4, -1]";
}

// The hierarchical rate prior integrates to 1 over the rates' range whatever the mean and shape, so the walk's bounds put
// no weight on them.
TEST(SampledLeakTest, TheTruncatedRatePriorIsADensityOnTheRatesRange) {
    for (const double mean : {0.005, 0.02, 0.2, 0.45}) {
        for (const double shape : {0.5, 2.0, 10.0, 150.0}) {
            // Midpoint rule in log x, fine enough for the 1/x-like behaviour near ERROR_RATE_LOWER at small shapes.
            constexpr int steps = 400000;
            const double a = std::log(ERROR_RATE_LOWER), b = std::log(ERROR_RATE_UPPER);
            double total = 0;
            for (int k = 0; k < steps; ++k) {
                const double x = std::exp(a + (k + 0.5) * (b - a) / steps);
                total += std::exp(truncatedGammaLogDensity(x, mean, shape)) * x * (b - a) / steps;
            }
            EXPECT_NEAR(total, 1.0, 1e-4) << "mean " << mean << ", shape " << shape;
        }
    }
    EXPECT_EQ(truncatedGammaLogDensity(0.6, 0.02, 2.0), -std::numeric_limits<double>::infinity());
}

// A restored leak from a Run with another prior is moved inside this one's range; one inside stays.
TEST(SampledLeakTest, ARestoredLeakOutsideThePriorIsMovedInside) {
    auto rng         = std::make_shared<boost::random::mt19937>(7);
    const auto study = leakStudy(*rng, 4);
    State state(study, rng, ModelChoice{.transmissionLeak = 0.01, .leakPrior = "loguniform:1e-4,0.1"});
    state.leakSize->initializeValue(0.3);
    state.keepLeakInsidePrior();
    EXPECT_LT(state.leakSize->value(), 0.1);
    EXPECT_GT(state.leakSize->value(), 0.099);
    state.leakSize->initializeValue(0.02);
    state.keepLeakInsidePrior();
    EXPECT_EQ(state.leakSize->value(), 0.02);
    EXPECT_EQ(canonicalLeakPrior("loguniform"), canonicalLeakPrior("loguniform:1e-7,0.1"));
    EXPECT_EQ(canonicalLeakPrior("exponential"), "exponential:100");
    EXPECT_EQ(canonicalLeakPrior("sideways"), "sideways");
    EXPECT_EQ(canonicalErrorModel("infection"), "");
    EXPECT_EQ(canonicalErrorModel(""), "");
    EXPECT_EQ(canonicalErrorModel("shared"), "shared");
    EXPECT_EQ(canonicalErrorModel("pooled"), "pooled");
}
