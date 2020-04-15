//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <Eigen/Core>

#include "core/parameters/Parameter.h"
#include "core/samplers/RandomWalkMH.h"

TEST(RandomWalkMHTest, NormalTest) {
    constexpr double TEST_MEAN = 5;
    constexpr double TEST_VARIANCE = 1;
    constexpr int TOTAL_DATA_POINTS = 100;

    struct NormalTestTarget {

        explicit NormalTestTarget(Parameter<double> &mean) : mean_(mean) {
            boost::random::normal_distribution<> dist{0, 1};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = dist(r) * TEST_VARIANCE + TEST_MEAN;
            }
        };

        double value() {
            boost::math::normal d(mean_.value(), 1);
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(boost::math::pdf(d, data_[i]));
            }
            return llik;
        }

        boost::random::mt19937 r;
        Parameter<double> &mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
    };


    Parameter<double> myMean(30);
    NormalTestTarget myTestTar(myMean);
    boost::random::mt19937 r;


    RandomWalkMH sampler(myMean, myTestTar, &r, 1);

    int i = 2000;
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
        results(i) = myMean.value();
    }

    auto resultsMean = results.mean();
    auto resultsStdDev = std::sqrt( (results - results.mean()).square().sum() / results.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Variance: " << sampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}

TEST(RandomWalkMHTest, DoubleWellTest) {
    struct DoubleWellTestTarget {

        explicit DoubleWellTestTarget(Parameter<double> &x, Parameter<double> &y) : x_(x), y_(y) {};

        double value() {
            return -(
                    (.25 * a_ * std::pow(x_.value(), 4)) -
                    (.5 * b_ * std::pow(x_.value(), 2)) +
                    (c_ * x_.value()) +
                    (.5 * d_ * std::pow(y_.value(), 2))
                    );
        }

        Parameter<double> &x_;
        Parameter<double> &y_;
        double a_ = 1;
        double b_ = 6;
        double c_ = 1;
        double d_ = 1;
    };

    Parameter<double> x(0);
    Parameter<double> y(0);
    DoubleWellTestTarget myTestTar(x, y);
    boost::random::mt19937 r;

    RandomWalkMH xSampler(x, myTestTar, &r, 10, 3, 100);
    RandomWalkMH ySampler(y, myTestTar, &r, 10, 3, 100);

    int i = 20000;
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
    std::cout << "Variance: " << xSampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << xSampler.acceptanceRate() << std::endl;
//    std::cout << xResults << std::endl;
//    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
//    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}