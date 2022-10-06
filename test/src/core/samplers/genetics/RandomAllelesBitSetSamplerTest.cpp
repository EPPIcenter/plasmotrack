//
// Created by Maxwell Murphy on 4/16/20.
//
//
#include <Eigen/Core>
#include <boost/random.hpp>
#include <utility>

#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::samplers;

TEST(RandomAllelesBitSetSamplerTest, AllelesBitSetTest) {
    using Alleles = AllelesBitSet<24>;

    struct AllelesBitSetTestTarget {
        explicit AllelesBitSetTestTarget(std::shared_ptr<Parameter<Alleles>> alleles) : alleles_(std::move(alleles)) {
            alleles_->add_post_change_listener([=, this]() {
                dirty = true;
            });
        }

        double value() {
            if (dirty) {
                value_ = 0;
                for (unsigned int i = 0; i < target.totalAlleles(); ++i) {
                    if (target.allele(i) == alleles_->value().allele(i)) {
                        value_ += 100;
                    }
                }
                dirty = false;
            }
            return value_;
        };

        [[nodiscard]] bool isDirty() const {
            return dirty;
        }

        Alleles target{"011010"};
        std::shared_ptr<Parameter<Alleles>> alleles_;
        bool dirty{true};
        double value_{0};
    };


    auto myAlleles = std::make_shared<Parameter<Alleles>>("000000");
    auto myTestTar = std::make_shared<AllelesBitSetTestTarget>(myAlleles);
    auto r         = std::make_shared<boost::random::mt19937>();

    genetics::RandomAllelesBitSetSampler sampler(myAlleles, myTestTar, r, 6);

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
        for (unsigned int j = 0; j < myAlleles->value().totalAlleles(); ++j) {
            results(j) += (myAlleles->value().allele(j)) / 1000.0;
        }
    }

    EXPECT_NEAR(results(0), 0, 1e5);
    EXPECT_NEAR(results(1), 1, 1e5);
    EXPECT_NEAR(results(2), 1, 1e5);
    EXPECT_NEAR(results(3), 0, 1e5);
    EXPECT_NEAR(results(4), 1, 1e5);
    EXPECT_NEAR(results(5), 0, 1e5);
}