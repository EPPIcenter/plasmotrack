//
// Created by Maxwell Murphy on 4/9/20.
//


#include <boost/random.hpp>
#include <Eigen/Core>

#include "gtest/gtest.h"

#include "core/datatypes/Simplex.h"
#include "core/io/serialize.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/general/SALTSampler.h"


using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::samplers;
using namespace transmission_nets::core::io;

TEST(SALTSamplerTest, SimplexTest) {

    struct SimplexTestTarget {
        explicit SimplexTestTarget(Parameter<Simplex> &freqs) : freqs(freqs) {
            freqs.add_post_change_listener([=, this]() {
                this->is_dirty = true;
            });
        }

        double value() {
            double llik = 0;
            for (unsigned int i = 0; i < target.size(); i++) {
                llik += target.at(i) * log(freqs.value().frequencies(i));
            }
            is_dirty = false;
            return llik;
        }

        [[nodiscard]] bool isDirty() const {
            return is_dirty;
        }

        std::vector<int> target{1, 2000, 3000, 4000};
        Parameter<Simplex> &freqs;
        bool is_dirty{true};
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
    int total_samples = 100000;
    i = total_samples;
    while (i > 0) {
        i--;
        sampler.update();
        for (unsigned int j = 0; j < mySimplex.value().totalElements(); ++j) {
            results(j) += mySimplex.value().frequencies(j) / total_samples;
        }
    }

    std::cout << serialize(std::vector<double>{results(0), results(1), results(2), results(3)}) << std::endl;
    EXPECT_NEAR(results(0), 1.0 / 9001, .015);
    EXPECT_NEAR(results(1), 2000.0 / 9001, .015);
    EXPECT_NEAR(results(2), 3000.0 / 9001, .015);
    EXPECT_NEAR(results(3), 4000.0 / 9001, .015);

}