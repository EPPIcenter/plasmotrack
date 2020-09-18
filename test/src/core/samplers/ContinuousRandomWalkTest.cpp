//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <Eigen/Core>

#include "core/parameters/Parameter.h"
#include "core/samplers/ContinuousRandomWalk.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::samplers;

TEST(ContinuousRandomWalkTest, NormalTest) {
    constexpr double TEST_MEAN = 5;
    constexpr double TEST_VARIANCE = 1;
    constexpr int TOTAL_DATA_POINTS = 100;

    struct NormalTestTarget {

        explicit NormalTestTarget(Parameter<double> &mean) : mean_(mean) {
            boost::random::normal_distribution<> dist{0, 1};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = dist(r) * TEST_VARIANCE + TEST_MEAN;
            }

            mean_.add_post_change_listener([=, this]() {
                dirty = true;
            });
        };

        bool isDirty() {
            return dirty;
        }

        double value() {
            if (dirty) {
                boost::math::normal d(mean_.value(), TEST_VARIANCE);
                value_ = 0;
                for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                    value_ += log(boost::math::pdf(d, data_[i]));
                }
                dirty = false;
                std::cout << value_ << " " << mean_.value() << std::endl;
            }

            return value_;
        }

        boost::random::mt19937 r;
        Parameter<double> &mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
        double value_;
        bool dirty{true};
    };


    Parameter<double> myMean(30);
    NormalTestTarget myTestTar(myMean);
    boost::random::mt19937 r;


    ContinuousRandomWalk sampler(myMean, myTestTar, &r, 3);
    sampler.setAdaptationRate(2);

    int i = 20000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
        std::cout << "Variance:" << sampler.variance() << std::endl;
    }

    Eigen::Array<double, 5000, 1> results;
    i = 5000;
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

TEST(ContinuousRandomWalkTest, DoubleWellTest) {
    struct DoubleWellTestTarget {

        DoubleWellTestTarget(Parameter<double> &x, Parameter<double> &y) : x_(x), y_(y) {
            x_.add_post_change_listener([=, this]() {
                this->dirty = true;
            });

            y_.add_post_change_listener([=, this]() {
                this->dirty = true;
            });
        };

        double value() {
            if (dirty) {
                value_ =  -(
                        (.25 * a_ * std::pow(x_.value(), 4)) -
                        (.5 * b_ * std::pow(x_.value(), 2)) +
                        (c_ * x_.value()) +
                        (.5 * d_ * std::pow(y_.value(), 2))
                );
                dirty = false;
            }
            return value_;
        }

        bool isDirty() {
            return dirty;
        }

        Parameter<double> &x_;
        Parameter<double> &y_;
        double a_ = 1;
        double b_ = 6;
        double c_ = 1;
        double d_ = 1;
        bool dirty = true;
        double value_;
    };

    Parameter<double> x(0);
    Parameter<double> y(0);
    DoubleWellTestTarget myTestTar(x, y);
    boost::random::mt19937 r;

    ContinuousRandomWalk xSampler(x, myTestTar, &r, 10, 3, 100);
    ContinuousRandomWalk ySampler(y, myTestTar, &r, 10, 3, 100);

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
}