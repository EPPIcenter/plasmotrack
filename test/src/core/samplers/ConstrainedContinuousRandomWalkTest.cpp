//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include "Eigen/Core"
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include <cmath>
#include <functional>
#include <utility>

#include "core/parameters/Parameter.h"
#include "core/abstract/observables/Dependency.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"

using namespace transmission_nets;
using namespace transmission_nets;

constexpr double TEST_PROB      = .15;
constexpr int TOTAL_DATA_POINTS = 100;
TEST(ConstrainedRandomWalkMHTest, BernoulliTest) {

    struct BernoulliTestTarget : public core::computation::Computation<core::computation::Likelihood>,
                                 public core::abstract::Observable<BernoulliTestTarget>,
                                 public core::abstract::Cacheable<BernoulliTestTarget>,
                                 public core::abstract::Checkpointable<BernoulliTestTarget, core::computation::Likelihood> {

        explicit BernoulliTestTarget(std::shared_ptr<core::parameters::Parameter<double>> prob) : prob_(std::move(prob)) {
            boost::random::bernoulli_distribution<> dist{TEST_PROB};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_.push_back(dist(r));
            }

            core::abstract::dependsOn(this, prob_, [this]() { dirty = true; });
            value_ = calculateValue(data_, prob_->value());
        };

        [[nodiscard]] bool isDirty() const {
            return dirty;
        }

        core::computation::Likelihood value() override {
            if (dirty) {
                value_ = calculateValue(data_, prob_->value());
                dirty  = false;
            }
            return value_;
        }

        [[nodiscard]] static core::computation::Likelihood calculateValue(const std::vector<double>& data, double prob) {
            boost::math::bernoulli d(prob);
            core::computation::Likelihood val = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                val += log(boost::math::pdf(d, data[i]));
            }
            return val;
        }

        boost::random::mt19937 r;
        std::shared_ptr<core::parameters::Parameter<double>> prob_;
        std::vector<double> data_;
        core::computation::Likelihood value_;
        bool dirty{true};
    };

    auto myProb    = std::make_shared<core::parameters::Parameter<double>>(.5);
    auto myTestTar = std::make_shared<BernoulliTestTarget>(myProb);
    auto r         = std::make_shared<boost::random::mt19937>();

    //    ConstrainedContinuousRandomWalk<0, 1> sampler2(myProb, myTestTar, &r, .01);
    core::samplers::ConstrainedContinuousRandomWalk<BernoulliTestTarget, boost::random::mt19937> sampler(myProb, myTestTar, 0, 1, r, .01);

    int i = 20000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    Eigen::Array<double, 1000, 1> results;
    i = 1000;
    while (i > 0) {
        i--;
        sampler.update();
        results(i) = myProb->value();
    }

    auto resultsMean   = results.mean();
    auto resultsStdDev = std::sqrt((results - results.mean()).square().sum() / results.size());
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Variance: " << sampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_PROB);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_PROB);
}
namespace {
    // A log density of the parameter, for checking the walk's stationary moments against known ones.
    struct DensityTarget {
        std::shared_ptr<core::parameters::Parameter<double>> x;
        std::function<double(double)> logDensity;

        double value() const { return logDensity(x->value()); }
        [[nodiscard]] bool isDirty() const { return false; }
    };

    std::pair<double, double> moments(const double lower, const double upper, const double start, std::function<double(double)> logDensity) {
        auto x      = std::make_shared<core::parameters::Parameter<double>>(start);
        auto target = std::make_shared<DensityTarget>(DensityTarget{x, std::move(logDensity)});
        auto rng    = std::make_shared<boost::random::mt19937>(20261008);
        core::samplers::ConstrainedContinuousRandomWalk<DensityTarget, boost::random::mt19937> sampler(x, target, lower, upper, rng, 1.0);
        for (int i = 0; i < 10000; ++i) {
            sampler.update();
        }
        constexpr int draws = 400000;
        double sum = 0, sumSquares = 0;
        for (int i = 0; i < draws; ++i) {
            sampler.update();
            sum += x->value();
            sumSquares += x->value() * x->value();
        }
        const double mean = sum / draws;
        return {mean, sumSquares / draws - mean * mean};
    }
}// namespace

// Weak targets, so that the Metropolis-Hastings correction for the logit transform matters.
TEST(ConstrainedRandomWalkMHTest, BetaMoments) {
    const auto [mean, variance] = moments(0, 1, 0.5, [](const double p) { return std::log(p) + 4 * std::log(1 - p); });
    EXPECT_NEAR(mean, 2.0 / 7, 0.005);
    EXPECT_NEAR(variance, 10.0 / (49 * 8), 0.001);
}

TEST(ConstrainedRandomWalkMHTest, FlatTargetOnWiderRange) {
    const auto [mean, variance] = moments(1, 20, 10, [](double) { return 0.0; });
    EXPECT_NEAR(mean, 10.5, 0.1);
    EXPECT_NEAR(variance, 19.0 * 19.0 / 12, 0.5);
}
