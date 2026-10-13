//
// SimplexLogRatioSampler must sample its target. The targets are Dirichlet densities on the simplex, whose means and
// variances are known; the sampler adapts first (scale and covariance), then runs fixed.
//

#include "gtest/gtest.h"

#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/general/SimplexLogRatioSampler.h"

#include <boost/random.hpp>

#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

using namespace transmission_nets::core;
using transmission_nets::core::samplers::general::SimplexLogRatioSampler;

namespace {
    struct DirichletTarget {
        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> freqs;
        std::vector<double> alpha;
        [[nodiscard]] double value() const {
            double out = 0;
            for (size_t i = 0; i < alpha.size(); ++i) {
                out += (alpha[i] - 1) * std::log(freqs->value().frequencies(i));
            }
            return out;
        }
    };
}// namespace

class SimplexLogRatioSamplerTest : public testing::TestWithParam<std::vector<double>> {};

INSTANTIATE_TEST_SUITE_P(Dirichlet, SimplexLogRatioSamplerTest,
                         testing::Values(std::vector<double>{2, 5, 10, 3}, std::vector<double>{1.5, 1, 4, 20, 2, 8, 3, 1.2, 6, 2.5}));

TEST_P(SimplexLogRatioSamplerTest, MatchesDirichletMoments) {
    const auto alpha = GetParam();
    const size_t k   = alpha.size();
    auto freqs       = std::make_shared<parameters::Parameter<datatypes::Simplex>>(std::vector<double>(k, 1.0 / k));
    auto target      = std::make_shared<DirichletTarget>(DirichletTarget{freqs, alpha});
    auto rng         = std::make_shared<boost::random::mt19937>(20261013);
    SimplexLogRatioSampler<DirichletTarget, boost::random::mt19937> sampler(freqs, target, rng, 1e-12);

    for (int i = 0; i < 20000; ++i) {
        sampler.update();
        sampler.adapt();
    }
    const double tuned = sampler.acceptanceRate();
    EXPECT_GT(tuned, 0.1);
    EXPECT_LT(tuned, 0.5);

    constexpr int draws = 400000;
    std::vector<double> mean(k, 0.0), meanSq(k, 0.0);
    for (int i = 0; i < draws; ++i) {
        sampler.update();
        for (size_t j = 0; j < k; ++j) {
            const double p = freqs->value().frequencies(j);
            mean[j] += p / draws;
            meanSq[j] += p * p / draws;
        }
    }
    const double a0 = std::accumulate(alpha.begin(), alpha.end(), 0.0);
    for (size_t j = 0; j < k; ++j) {
        const double m = alpha[j] / a0;
        const double v = alpha[j] * (a0 - alpha[j]) / (a0 * a0 * (a0 + 1));
        EXPECT_NEAR(mean[j], m, 0.05 * std::sqrt(v) + 0.002) << "frequency " << j;
        EXPECT_NEAR(meanSq[j] - mean[j] * mean[j], v, 0.05 * v + 1e-5) << "frequency " << j;
    }
}
