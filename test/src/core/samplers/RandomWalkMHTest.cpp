//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"
#include "Eigen/Core"

#include "core/parameters/Parameter.h"
#include "core/samplers/RandomWalkMH.h"

constexpr double TEST_MEAN = 5;
constexpr int TOTAL_DATA_POINTS = 100;
TEST(RandomWalkMHTest, NormalTest) {

    struct NormalTestTarget {

        explicit NormalTestTarget(Parameter<double> &mean) : mean_(mean) {
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = gsl_ran_gaussian(r, 3) + TEST_MEAN;
            }
        };

        double value() {
            double llik = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                llik += log(gsl_ran_gaussian_pdf(mean_.value() - data_(i), 1));
            }
            return llik;
        }

        gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        Parameter<double> &mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
    };


    Parameter<double> myMean(30);
    NormalTestTarget myTestTar(myMean);
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, 400);

    RandomWalkMH sampler(myMean, myTestTar, r, 1);

    int i = 5000;
    while (i > 0) {
        i--;
        sampler.update();
//        std::cout << sampler.variance() << std::endl;
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