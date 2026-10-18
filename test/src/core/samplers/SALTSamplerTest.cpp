//
// Created by Maxwell Murphy on 4/9/20.
//
//
//
#include <Eigen/Core>
#include <boost/random.hpp>

#include <array>
#include <cmath>

#include "gtest/gtest.h"

#include "core/datatypes/Simplex.h"
#include "core/io/serialize.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/general/SALTSampler.h"


using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::samplers;
using namespace transmission_nets::core::io;

TEST(SALTSamplerTest, SimplexTest) {

    struct SimplexTestTarget {
        explicit SimplexTestTarget(std::shared_ptr<Parameter<Simplex>> freqs) : freqs_(std::move(freqs)) {
            freqs_->add_post_change_listener([=, this]() {
                this->is_dirty = true;
            });
        }

        double value() {
            double llik = 0;
            for (unsigned int i = 0; i < target.size(); i++) {
                llik += target.at(i) * log(freqs_->value().frequencies(i));
            }
            is_dirty = false;
            return llik;
        }

        [[nodiscard]] bool isDirty() const {
            return is_dirty;
        }

        std::vector<int> target{1, 2000, 3000, 4000};
        std::shared_ptr<Parameter<Simplex>> freqs_;
        bool is_dirty{true};
    };

    auto mySimplex = std::make_shared<Parameter<Simplex>>(std::vector{.1, .1, .1, .1});
    auto st        = std::make_shared<SimplexTestTarget>(mySimplex);
    auto r         = std::make_shared<boost::random::mt19937>();

    SALTSampler sampler(mySimplex, st, r, 1, .1, 10);

    int i = 50000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    Eigen::Array<double, 4, 1> results;
    results.setZero();
    int total_samples = 100000;
    i                 = total_samples;
    while (i > 0) {
        i--;
        sampler.update();
        for (unsigned int j = 0; j < mySimplex->value().totalElements(); ++j) {
            results(j) += mySimplex->value().frequencies(j) / total_samples;
        }
    }

    std::cout << serialize(std::vector<double>{results(0), results(1), results(2), results(3)}) << std::endl;
    std::cout << "Acceptances: ";
    for (int j = 0; j < 4; ++j) {
        std::cout << sampler.acceptanceRate(j) << ", ";
    }
    std::cout << std::endl;
    std::cout << "Variances: ";
    for (int j = 0; j < 4; ++j) {
        std::cout << sampler.getVariance(j) << ", ";
    }
    std::cout << std::endl;
    EXPECT_NEAR(results(0), 1.0 / 9001, .015);
    EXPECT_NEAR(results(1), 2000.0 / 9001, .015);
    EXPECT_NEAR(results(2), 3000.0 / 9001, .015);
    EXPECT_NEAR(results(3), 4000.0 / 9001, .015);
}
// A weak target, so that the Metropolis-Hastings correction matters: under a Dirichlet(2, 3, 5) the frequencies have
// means alpha / 10 and variances alpha (10 - alpha) / 1100.
TEST(SALTSamplerTest, DirichletMoments) {
    struct DirichletTarget {
        explicit DirichletTarget(std::shared_ptr<Parameter<Simplex>> freqs) : freqs_(std::move(freqs)) {
            freqs_->add_post_change_listener([this]() { dirty_ = true; });
        }

        double value() {
            double out = 0;
            for (unsigned int i = 0; i < alpha.size(); i++) {
                out += (alpha[i] - 1) * std::log(freqs_->value().frequencies(i));
            }
            dirty_ = false;
            return out;
        }

        [[nodiscard]] bool isDirty() const { return dirty_; }

        std::vector<double> alpha{2, 3, 5};
        std::shared_ptr<Parameter<Simplex>> freqs_;
        bool dirty_ = true;
    };

    auto freqs  = std::make_shared<Parameter<Simplex>>(std::vector{1.0 / 3, 1.0 / 3, 1.0 / 3});
    auto target = std::make_shared<DirichletTarget>(freqs);
    auto rng    = std::make_shared<boost::random::mt19937>(20261008);

    SALTSampler sampler(freqs, target, rng, 1);
    sampler.setLowerLimit_(1e-9);

    for (int i = 0; i < 10000; ++i) {
        sampler.update();
    }
    constexpr int draws = 400000;
    std::array<double, 3> sum{};
    std::array<double, 3> sumSquares{};
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (unsigned int j = 0; j < 3; ++j) {
            const double f = freqs->value().frequencies(j);
            sum[j] += f;
            sumSquares[j] += f * f;
        }
    }
    const std::vector<double>& alpha = target->alpha;
    for (unsigned int j = 0; j < 3; ++j) {
        const double mean     = sum[j] / draws;
        const double variance = sumSquares[j] / draws - mean * mean;
        EXPECT_NEAR(mean, alpha[j] / 10, 0.005) << "frequency " << j;
        EXPECT_NEAR(variance, alpha[j] * (10 - alpha[j]) / 1100, 0.001) << "frequency " << j;
    }
}

namespace {
    struct WeakDirichletTarget {
        explicit WeakDirichletTarget(std::shared_ptr<Parameter<Simplex>> freqs, std::vector<double> alpha = {2, 3, 5}) : alpha(std::move(alpha)), freqs_(std::move(freqs)) {
            freqs_->add_post_change_listener([this]() { dirty_ = true; });
        }

        double value() {
            double out = 0;
            for (unsigned int i = 0; i < alpha.size(); i++) {
                out += (alpha[i] - 1) * std::log(freqs_->value().frequencies(i));
            }
            dirty_ = false;
            return out;
        }

        [[nodiscard]] bool isDirty() const { return dirty_; }

        std::vector<double> alpha;
        std::shared_ptr<Parameter<Simplex>> freqs_;
        bool dirty_ = true;
    };
}// namespace

// With the Dirichlet move, fitted and then fixed, the Dirichlet(2, 3, 5) moments still hold, and the fitted move is
// accepted often.
TEST(SALTSamplerTest, DirichletMoveKeepsDirichletMoments) {
    auto freqs  = std::make_shared<Parameter<Simplex>>(std::vector{1.0 / 3, 1.0 / 3, 1.0 / 3});
    auto target = std::make_shared<WeakDirichletTarget>(freqs);
    auto rng    = std::make_shared<boost::random::mt19937>(20261017);

    SALTSampler sampler(freqs, target, rng, 1);
    sampler.setLowerLimit_(1e-9);
    sampler.enableDirichletMove({200, 1000}, 0.25, 10);

    for (int i = 0; i < 5000; ++i) {
        sampler.update();
    }
    ASSERT_EQ(sampler.dirichletAlpha().size(), 3u);
    const auto fittedCounts = sampler.dirichletCounts();
    constexpr int draws     = 400000;
    std::array<double, 3> sum{};
    std::array<double, 3> sumSquares{};
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (unsigned int j = 0; j < 3; ++j) {
            const double f = freqs->value().frequencies(j);
            sum[j] += f;
            sumSquares[j] += f * f;
        }
    }
    const std::vector<double>& alpha = target->alpha;
    for (unsigned int j = 0; j < 3; ++j) {
        const double mean     = sum[j] / draws;
        const double variance = sumSquares[j] / draws - mean * mean;
        EXPECT_NEAR(mean, alpha[j] / 10, 0.005) << "frequency " << j;
        EXPECT_NEAR(variance, alpha[j] * (10 - alpha[j]) / 1100, 0.001) << "frequency " << j;
    }
    const auto counts = sampler.dirichletCounts();
    EXPECT_GT(static_cast<double>(counts.accepted - fittedCounts.accepted), 0.5 * static_cast<double>(counts.proposed - fittedCounts.proposed));
}

// Moved in blocks of 3 (3, 3 and 4 of the 10 frequencies, partitioned afresh each move), a Dirichlet target's moments
// still hold.
TEST(SALTSamplerTest, DirichletBlockMoveKeepsDirichletMoments) {
    const std::vector<double> alpha{2, 3, 5, 1.5, 4, 2.5, 6, 3.5, 2, 4.5};
    const double total = 34;
    auto freqs         = std::make_shared<Parameter<Simplex>>(std::vector<double>(alpha.size(), 1.0 / alpha.size()));
    auto target        = std::make_shared<WeakDirichletTarget>(freqs, alpha);
    auto rng           = std::make_shared<boost::random::mt19937>(20261018);

    SALTSampler sampler(freqs, target, rng, 1);
    sampler.setLowerLimit_(1e-9);
    sampler.enableDirichletMove({200, 1000}, 0.25, 10, 1, 3);

    for (int i = 0; i < 5000; ++i) {
        sampler.update();
    }
    ASSERT_EQ(sampler.dirichletAlpha().size(), alpha.size());
    const auto fittedCounts = sampler.dirichletCounts();
    constexpr int draws     = 200000;
    std::vector<double> sum(alpha.size()), sumSquares(alpha.size());
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (unsigned int j = 0; j < alpha.size(); ++j) {
            const double f = freqs->value().frequencies(j);
            sum[j] += f;
            sumSquares[j] += f * f;
        }
    }
    for (unsigned int j = 0; j < alpha.size(); ++j) {
        const double mean     = sum[j] / draws;
        const double variance = sumSquares[j] / draws - mean * mean;
        EXPECT_NEAR(mean, alpha[j] / total, 0.003) << "frequency " << j;
        EXPECT_NEAR(variance, alpha[j] * (total - alpha[j]) / (total * total * (total + 1)), 0.0003) << "frequency " << j;
    }
    const auto counts = sampler.dirichletCounts();
    EXPECT_EQ(counts.proposed - fittedCounts.proposed, 3u * draws);
    EXPECT_GT(static_cast<double>(counts.accepted - fittedCounts.accepted), 0.5 * static_cast<double>(counts.proposed - fittedCounts.proposed));
}

// A sampler restored mid-window continues exactly as the original does; a state without a fit restores and fits
// afresh, and a fit of the wrong size is refused.
TEST(SALTSamplerTest, DirichletFitRestores) {
    auto make = [](std::shared_ptr<Parameter<Simplex>> freqs, std::shared_ptr<boost::random::mt19937> rng) {
        auto target  = std::make_shared<WeakDirichletTarget>(freqs);
        auto sampler = std::make_unique<SALTSampler<WeakDirichletTarget>>(freqs, target, rng, 1);
        sampler->setLowerLimit_(1e-9);
        sampler->enableDirichletMove({20, 50}, 0.25, 10);
        return std::pair{std::move(sampler), target};
    };
    auto freqsA       = std::make_shared<Parameter<Simplex>>(std::vector{1.0 / 3, 1.0 / 3, 1.0 / 3});
    auto rngA         = std::make_shared<boost::random::mt19937>(7);
    auto [a, targetA] = make(freqsA, rngA);
    for (int i = 0; i < 35; ++i) {
        a->update();
    }
    const auto saved = a->adaptationState();
    ASSERT_FALSE(saved.fit.empty());

    auto freqsB       = std::make_shared<Parameter<Simplex>>(freqsA->value().frequencies());
    auto rngB         = std::make_shared<boost::random::mt19937>(*rngA);
    auto [b, targetB] = make(freqsB, rngB);
    b->restoreAdaptationState(saved);
    for (int i = 0; i < 100; ++i) {
        a->update();
        b->update();
    }
    EXPECT_EQ(freqsB->value().frequencies(), freqsA->value().frequencies());
    EXPECT_EQ(b->dirichletAlpha(), a->dirichletAlpha());
    EXPECT_EQ(b->adaptationState(), a->adaptationState());

    auto [c, targetC] = make(std::make_shared<Parameter<Simplex>>(std::vector{1.0 / 3, 1.0 / 3, 1.0 / 3}), std::make_shared<boost::random::mt19937>(1));
    auto withoutFit   = saved;
    withoutFit.fit.clear();
    c->restoreAdaptationState(withoutFit);
    EXPECT_TRUE(c->dirichletAlpha().empty());
    auto wrongSize = saved;
    wrongSize.fit.pop_back();
    EXPECT_THROW(c->restoreAdaptationState(wrongSize), std::invalid_argument);
}
