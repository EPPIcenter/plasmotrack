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
