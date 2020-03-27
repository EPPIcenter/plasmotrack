//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"
#include "Eigen/Core"

#include "core/parameters/Parameter.h"
#include "core/samplers/ConstrainedRandomWalkMH.h"

constexpr double TEST_PROB = .15;
constexpr int TOTAL_DATA_POINTS = 100;
TEST(ConstrainedRandomWalkMHTest, BernoulliTest) {

    struct BernoulliTestTarget {

        explicit BernoulliTestTarget(Parameter<double> &prob) : prob_(prob) {
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_.push_back(gsl_ran_bernoulli(r, TEST_PROB));
            }
        };

        double value() {
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(gsl_ran_bernoulli_pdf(data_[i], prob_.value()));
            }
            return llik;
        }

    private:
        gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        Parameter<double> &prob_;
        std::vector<double> data_;
    };

    Parameter<double> myProb(.5);
    BernoulliTestTarget myTestTar(myProb);
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);

    ConstrainedRandomWalkMH sampler(myProb, myTestTar, r, .01, 0, 1);

    int i = 10000;
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