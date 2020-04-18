//
// Created by Maxwell Murphy on 4/9/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "Eigen/Core"

#include "core/parameters/Parameter.h"

#include "core/datatypes/Simplex.h"

#include "core/samplers/SALTSampler.h"


TEST(SALTSamplerTest, SimplexTest) {

    struct SimplexTestTarget {
        SimplexTestTarget(Parameter<Simplex> &freqs) : freqs(freqs) {}

        double value() {
            double llik = 0;
            for (unsigned int i = 0; i < target.size(); i++) {
                llik += target.at(i) * log(freqs.value().frequencies(i));
            }
            return llik;
        }

        std::vector<int> target{1, 2000, 3000, 4000};
        Parameter<Simplex> &freqs;
    };

    Parameter<Simplex> mySimplex{.1, .1, .1, .1};
    SimplexTestTarget st(mySimplex);
    boost::random::mt19937 r;

    SALTSampler sampler(mySimplex, st, &r);

    int i = 5000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    Eigen::Array<double, 4, 1> results;
    results.setZero();
    i = 1000;
    while (i > 0) {
        i--;
        sampler.update();
        for (unsigned int j = 0; j < mySimplex.value().totalElements(); ++j) {
            results(j) += mySimplex.value().frequencies(j) / 1000.0;
        }
    }

    EXPECT_NEAR(results(0), 1.0 / 9001, .015);
    EXPECT_NEAR(results(1), 2000.0 / 9001, .015);
    EXPECT_NEAR(results(2), 3000.0 / 9001, .015);
    EXPECT_NEAR(results(3), 4000.0 / 9001, .015);

}