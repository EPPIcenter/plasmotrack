//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "Eigen/Core"

#include "core/parameters/Parameter.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::samplers;

constexpr double TEST_PROB = .15;
constexpr int TOTAL_DATA_POINTS = 100;
TEST(ConstrainedRandomWalkMHTest, BernoulliTest) {

    struct BernoulliTestTarget {

        explicit BernoulliTestTarget(Parameter<double> &prob) : prob_(prob) {
            boost::random::bernoulli_distribution<> dist{TEST_PROB};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_.push_back(dist(r));
            }
        };

        bool isDirty() {
            return true;
        }

        double value() {
            boost::math::bernoulli d(prob_.value());
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(boost::math::pdf(d, data_[i]));
            }
            return llik;
        }

        boost::random::mt19937 r;
        Parameter<double> &prob_;
        std::vector<double> data_;
    };

    Parameter<double> myProb(.5);
    BernoulliTestTarget myTestTar(myProb);
    boost::random::mt19937 r;

//    ConstrainedContinuousRandomWalk<0, 1> sampler2(myProb, myTestTar, &r, .01);
    ConstrainedContinuousRandomWalk<BernoulliTestTarget, boost::random::mt19937> sampler(myProb, myTestTar, 0, 1, &r, .01);

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
        results(i) = myProb.value();
    }

    auto resultsMean = results.mean();
    auto resultsStdDev = std::sqrt( (results - results.mean()).square().sum() / results.size() );
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Variance: " << sampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_PROB);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_PROB);
}