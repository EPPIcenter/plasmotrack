//
// Created by Maxwell Murphy on 4/16/20.
//

#include <boost/random.hpp>
#include <Eigen/Core>

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"
#include "core/datatypes/Alleles.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"

TEST(RandomAllelesBitSetSamplerTest, AllelesBitSetTest) {
    using Alleles = AllelesBitSet<24>;

    struct AllelesBitSetTestTarget {
        explicit AllelesBitSetTestTarget(Parameter<Alleles> &alleles) : alleles_(alleles) {}

        double value() {
            double llik = 0;
            for (unsigned int i = 0; i < target.totalAlleles(); ++i) {
                if(target.allele(i) == alleles_.value().allele(i)) {
                    llik += 10;
                }
            }
            return llik;
        };

        Alleles target{"011010"};
        Parameter<Alleles> &alleles_;

    };


    Parameter<Alleles> myAlleles("000000");
    AllelesBitSetTestTarget myTestTar(myAlleles);
    boost::random::mt19937 r;

    RandomAllelesBitSetSampler sampler(myAlleles, myTestTar, &r);

    int i = 5000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    Eigen::Array<double, 6, 1> results;
    i = 1000;
    while (i > 0) {
        i--;
        sampler.update();
        for (unsigned int j = 0; j < myAlleles.value().totalAlleles(); ++j) {
            results(j) += (myAlleles.value().allele(j)) / 1000.0;
        }
    }

    EXPECT_NEAR(results(0), 0, 1e5);
    EXPECT_NEAR(results(1), 1, 1e5);
    EXPECT_NEAR(results(2), 1, 1e5);
    EXPECT_NEAR(results(3), 0, 1e5);
    EXPECT_NEAR(results(4), 1, 1e5);
    EXPECT_NEAR(results(5), 0, 1e5);
}