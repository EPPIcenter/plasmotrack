//
// Created by Maxwell Murphy on 4/16/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <Eigen/Core>

#include "core/parameters/Parameter.h"
#include "core/samplers/DiscreteRandomWalk.h"

TEST(DiscreteRandomWalkTest, NormalTest) {
    constexpr double TEST_MEAN = 5.51;
    constexpr double TEST_VARIANCE = 10;
    constexpr int TOTAL_DATA_POINTS = 30;

    struct NormalTestTarget {

        explicit NormalTestTarget(Parameter<int> &mean) : mean_(mean) {
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = dist(r) * TEST_VARIANCE + TEST_MEAN;
            }
        };

        bool isDirty() {
            return true;
        }

        double value() {
            boost::math::normal d(mean_.value(), TEST_VARIANCE);
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(boost::math::pdf(d, data_[i]));
            }
            return llik;
        }

        boost::random::normal_distribution<> dist{0, 1};
        boost::random::mt19937 r;
        Parameter<int> &mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
    };


    Parameter<int> myMean(30);
    NormalTestTarget myTestTar(myMean);
    boost::random::mt19937 r;

    DiscreteRandomWalk sampler(myMean, myTestTar, &r, 3);

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
        results(i) = myMean.value();
    }

    auto resultsMean = results.mean();
    auto resultsStdDev = std::sqrt( (results - results.mean()).square().sum() / results.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}

TEST(DiscreteRandomWalkTest, DoubleWellTest) {
    struct DoubleWellTestTarget {

        explicit DoubleWellTestTarget(Parameter<int> &x, Parameter<int> &y) : x_(x), y_(y) {};

        double value() {
            return -(
                    (.25 * a_ * std::pow(x_.value(), 4)) -
                    (.5 * b_ * std::pow(x_.value(), 2)) +
                    (c_ * x_.value()) +
                    (.5 * d_ * std::pow(y_.value(), 2))
            );
        }

        Parameter<int> &x_;
        Parameter<int> &y_;
        double a_ = 1;
        double b_ = 6;
        double c_ = 1;
        double d_ = 1;
    };

    Parameter<int> x(0);
    Parameter<int> y(0);
    DoubleWellTestTarget myTestTar(x, y);
    boost::random::mt19937 r;

    DiscreteRandomWalk xSampler(x, myTestTar, &r, 3);
    DiscreteRandomWalk ySampler(y, myTestTar, &r, 3);

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
        xResults(i) = x.value();
        yResults(i) = y.value();
    }

    auto resultsMean = xResults.mean();
    auto resultsStdDev = std::sqrt( (xResults - xResults.mean()).square().sum() / xResults.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Acceptance Ratio: " << xSampler.acceptanceRate() << std::endl;
//    std::cout << xResults << std::endl;
//    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
//    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}