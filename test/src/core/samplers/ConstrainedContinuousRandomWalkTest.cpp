//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "Eigen/Core"

#include "core/parameters/Parameter.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"

using namespace transmission_nets;
using namespace transmission_nets;

constexpr double TEST_PROB = .15;
constexpr int TOTAL_DATA_POINTS = 100;
TEST(ConstrainedRandomWalkMHTest, BernoulliTest) {

    struct BernoulliTestTarget : public core::computation::Computation<core::computation::Likelihood>,
                                 public core::abstract::Observable<BernoulliTestTarget>,
                                 public core::abstract::Cacheable<BernoulliTestTarget>,
                                 public core::abstract::Checkpointable<BernoulliTestTarget, core::computation::Likelihood>  {

        explicit BernoulliTestTarget(std::shared_ptr<core::parameters::Parameter<double>> prob) : prob_(std::move(prob)) {
            boost::random::bernoulli_distribution<> dist{TEST_PROB};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_.push_back(dist(r));
            }

            prob_->registerCacheableCheckpointTarget(this);
            prob_->add_post_change_listener([=, this]() {
                dirty = true;
            });
            value_ = calculateValue(data_, prob_->value());
        };

        [[nodiscard]] bool isDirty() const {
            return dirty;
        }

        core::computation::Likelihood value() override {
            if (dirty) {
                value_ = calculateValue(data_, prob_->value());
                dirty = false;
            }
            return value_;
        }

        [[nodiscard]] static core::computation::Likelihood calculateValue(const std::vector<double> &data, double prob) {
            boost::math::bernoulli d(prob);
            core::computation::Likelihood val = 0;
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                val += log(boost::math::pdf(d, data[i]));
            }
            return val;
        }

        boost::random::mt19937 r;
        std::shared_ptr<core::parameters::Parameter<double>> prob_;
        std::vector<double> data_;
        core::computation::Likelihood value_;
        bool dirty{true};
    };

    auto myProb = std::make_shared<core::parameters::Parameter<double>>(.5);
    auto myTestTar = std::make_shared<BernoulliTestTarget>(myProb);
    auto r = std::make_shared<boost::random::mt19937>();

//    ConstrainedContinuousRandomWalk<0, 1> sampler2(myProb, myTestTar, &r, .01);
    core::samplers::ConstrainedContinuousRandomWalk<BernoulliTestTarget, boost::random::mt19937> sampler(myProb, myTestTar, 0, 1, r, .01);

    int i = 20000;
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
        results(i) = myProb->value();
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