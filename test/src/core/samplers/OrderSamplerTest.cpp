//
// Created by Maxwell Murphy on 4/14/20.
//

#include <boost/random.hpp>
#include <Eigen/Core>

#include "gtest/gtest.h"

#include "core/parameters/Ordering.h"
#include "core/samplers/topology/OrderSampler.h"


using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::samplers;
using namespace transmission_nets::core::samplers::topology;

TEST(OrderSamplerTest, OrderTest) {
    struct OrderingTestTarget {

        explicit OrderingTestTarget(Ordering<int> &ordering) : ordering_(ordering) {
            ordering_.add_post_change_listener([=, this]() {
                dirty = true;
            });
        };

        double value() {
            if (dirty) {
                value_ = 0;
                for (unsigned int j = 0; j < ordering_.value().size() - 1; ++j) {
                    for (unsigned int k = j + 1; k < ordering_.value().size(); ++k) {
                        if (*(ordering_.value().at(j)) < *(ordering_.value().at(k))) {
                            value_ += 10;
                        } else {
                            value_ -= 10;
                        }
                    }
                }
                dirty = false;
            }
            return value_;
        }

        bool isDirty() {
            return dirty;
        }

        Ordering<int> &ordering_;
        bool dirty{true};
        double value_{0};
    };

    auto a = std::make_shared<int>(1);
    auto b = std::make_shared<int>(2);
    auto c = std::make_shared<int>(3);
    auto d = std::make_shared<int>(4);
    auto e = std::make_shared<int>(5);
    auto f = std::make_shared<int>(6);
    auto g = std::make_shared<int>(7);


    Ordering<int> myOrdering({d, b, c, a, e, g, f});
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