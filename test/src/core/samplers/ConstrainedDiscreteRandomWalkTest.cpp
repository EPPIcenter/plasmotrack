//
// Created by Maxwell Murphy on 4/16/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <Eigen/Core>

#include "core/parameters/Parameter.h"
#include "core/samplers/general/ConstrainedDiscreteRandomWalk.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::samplers;

TEST(ConstrainedDiscreteRandomWalkTest, NormalTest) {
    constexpr double TEST_MEAN = 5.51;
    constexpr double TEST_VARIANCE = 10;
    constexpr int TOTAL_DATA_POINTS = 30;

    struct NormalTestTarget {

        explicit NormalTestTarget(std::shared_ptr<Parameter<int>> mean) : mean_(std::move(mean)) {
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = dist(r) * TEST_VARIANCE + TEST_MEAN;
            }
        };

        bool isDirty() {
            return true;
        }

        double value() {
            boost::math::normal d(mean_->value(), TEST_VARIANCE);
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(boost::math::pdf(d, data_[i]));
            }
            return llik;
        }

        boost::random::normal_distribution<> dist{0, 1};
        boost::random::mt19937 r;
        std::shared_ptr<Parameter<int>> mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
    };


    auto myMean = std::make_shared<Parameter<int>>(5);
    auto myTestTar = std::make_shared<NormalTestTarget>(myMean);
    auto r = std::make_shared<boost::random::mt19937>();

    ConstrainedDiscreteRandomWalk<3, 11, NormalTestTarget, boost::random::mt19937> sampler(myMean, myTestTar, r, 3);

    int i = 2000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    constexpr int total_samples = 2000;
    Eigen::Array<double, total_samples, 1> results;
    i = total_samples;
    while (i > 0) {
        i--;
        sampler.update();
        results(i) = myMean->value();
    }

    auto resultsMean = results.mean();
    auto resultsStdDev = std::sqrt( (results - results.mean()).square().sum() / results.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}

TEST(ConstrainedDiscreteRandomWalkTest, DoubleWellTest) {
    struct DoubleWellTestTarget {

        explicit DoubleWellTestTarget(std::shared_ptr<Parameter<int>> x, std::shared_ptr<Parameter<int>> y) : x_(std::move(x)), y_(std::move(y)) {};

        double value() {
            return -(
                    (.25 * a_ * std::pow(x_->value(), 4)) -
                    (.5 * b_ * std::pow(x_->value(), 2)) +
                    (c_ * x_->value()) +
                    (.5 * d_ * std::pow(y_->value(), 2))
            );
        }

        bool isDirty() {
            return true;
        }

        std::shared_ptr<Parameter<int>> x_;
        std::shared_ptr<Parameter<int>> y_;
        double a_ = 1;
        double b_ = 6;
        double c_ = 1;
        double d_ = 1;
    };

    auto x = std::make_shared<Parameter<int>>(0);
    auto y = std::make_shared<Parameter<int>>(0);
    auto myTestTar = std::make_shared<DoubleWellTestTarget>(x, y);
    auto r = std::make_shared<boost::random::mt19937>();

    ConstrainedDiscreteRandomWalk< -5, 5, DoubleWellTestTarget, boost::random::mt19937> xSampler(x, myTestTar, r, 3);
    ConstrainedDiscreteRandomWalk< -5, 5, DoubleWellTestTarget, boost::random::mt19937> ySampler(y, myTestTar, r, 3);

    int i = 2000;
    while (i > 0) {
        i--;
        xSampler.update();
        xSampler.adapt();
        ySampler.update();
        ySampler.adapt();
    }

    constexpr int total_samples = 2000;
    Eigen::Array<double, total_samples, 1> xResults;
    Eigen::Array<double, total_samples, 1> yResults;
    i = total_samples;
    while (i > 0) {
        i--;
        xSampler.update();
        ySampler.update();
        xResults(i) = x->value();
        yResults(i) = y->value();
    }

    auto resultsMean = xResults.mean();
    auto resultsStdDev = std::sqrt( (xResults - xResults.mean()).square().sum() / xResults.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Acceptance Ratio: " << xSampler.acceptanceRate() << std::endl;
}