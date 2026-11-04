//
// DirichletIndependenceSampler must leave its target invariant, with a proposal centred away from the target's mean.
//

#include "gtest/gtest.h"

#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/general/DirichletIndependenceSampler.h"

#include <boost/random.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::samplers;

namespace {
    // Under a Dirichlet(2, 3, 5) the frequencies have means alpha / 10 and variances alpha (10 - alpha) / 1100.
    struct DirichletTarget {
        explicit DirichletTarget(std::shared_ptr<Parameter<Simplex>> freqs) : freqs_(std::move(freqs)) {
            freqs_->add_post_change_listener([this]() { dirty_ = true; });
        }

        double value() {
            double out = 0;
            for (unsigned int i = 0; i < alpha.size(); i++) {
                out += (alpha[i] - 1) * std::log(freqs_->value().frequencies(i));
            }
            dirty_ = false;
            return out;
        }

        [[nodiscard]] bool isDirty() const { return dirty_; }

        std::vector<double> alpha{2, 3, 5};
        std::shared_ptr<Parameter<Simplex>> freqs_;
        bool dirty_ = true;
    };
}// namespace

TEST(DirichletIndependenceSamplerTest, LeavesTargetInvariant) {
    auto freqs  = std::make_shared<Parameter<Simplex>>(std::vector{1.0 / 3, 1.0 / 3, 1.0 / 3});
    auto target = std::make_shared<DirichletTarget>(freqs);
    auto rng    = std::make_shared<boost::random::mt19937>(20261017);

    // Centred on (0.5, 0.3, 0.2), not the target's (0.2, 0.3, 0.5).
    DirichletIndependenceSampler sampler(freqs, target, rng, {10, 6, 4}, 0.25, 1e-9);

    for (int i = 0; i < 10000; ++i) {
        sampler.update();
    }
    constexpr int draws = 400000;
    std::array<double, 3> sum{};
    std::array<double, 3> sumSquares{};
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (unsigned int j = 0; j < 3; ++j) {
            const double f = freqs->value().frequencies(j);
            sum[j] += f;
            sumSquares[j] += f * f;
        }
    }
    const std::vector<double>& alpha = target->alpha;
    for (unsigned int j = 0; j < 3; ++j) {
        const double mean     = sum[j] / draws;
        const double variance = sumSquares[j] / draws - mean * mean;
        EXPECT_NEAR(mean, alpha[j] / 10, 0.005) << "frequency " << j;
        EXPECT_NEAR(variance, alpha[j] * (10 - alpha[j]) / 1100, 0.001) << "frequency " << j;
    }
    EXPECT_GT(sampler.counts().accepted, sampler.counts().proposed / 20);
}
