//
// Created by Maxwell Murphy on 4/14/20.
//

#include <boost/random.hpp>
#include <Eigen/Core>

#include "gtest/gtest.h"

#include "core/parameters/Ordering.h"
#include "core/samplers/OrderSampler.h"



TEST(OrderSamplerTest, OrderTest) {
    struct OrderingTestTarget {

        explicit OrderingTestTarget(Ordering<int> &ordering) : ordering_(ordering) {};

        double value() {
            double llik = 0;
            for (unsigned int j = 0; j < ordering_.value().size() - 1; ++j) {
                for (unsigned int k = j + 1; k < ordering_.value().size(); ++k) {
                    if (*(ordering_.value().at(j)) < *(ordering_.value().at(k))) {
                        llik += 10;
                    } else {
                        llik -= 10;
                    }
                }
            }
            return llik;
        }

        bool isDirty() {
            return true;
        }

        Ordering<int> &ordering_;
    };

    int a = 1;
    int b = 2;
    int c = 3;
    int d = 4;
    int e = 5;
    int f = 6;
    int g = 7;


    Ordering<int> myOrdering({&d, &b, &c, &a, &e, &g, &f});
    OrderingTestTarget myTestTar(myOrdering);
    boost::random::mt19937 r;

    OrderSampler sampler(myOrdering, myTestTar, &r, 1);

    int i = 500;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
    }

    Eigen::Array<double, 7, 1> results;
    results.resize(myOrdering.value().size());
    results.setZero();
    i = 1000;
    while (i > 0) {
        i--;
        sampler.update();
        for (unsigned int j = 0; j < myOrdering.value().size(); ++j) {
            results(j) += *(myOrdering.value().at(j)) / 1000.0;
        }
    }

    EXPECT_NEAR(results(0), 1, 1e5);
    EXPECT_NEAR(results(1), 2, 1e5);
    EXPECT_NEAR(results(2), 3, 1e5);
    EXPECT_NEAR(results(3), 4, 1e5);
    EXPECT_NEAR(results(4), 5, 1e5);
    EXPECT_NEAR(results(5), 6, 1e5);
    EXPECT_NEAR(results(6), 7, 1e5);

}