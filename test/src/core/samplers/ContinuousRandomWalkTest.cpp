//
// Created by Maxwell Murphy on 3/25/20.
//

#include "gtest/gtest.h"

#include <Eigen/Core>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"
#include "core/samplers/general/ContinuousRandomWalk.h"

//using namespace transmission_nets::core::parameters;
//using namespace transmission_nets::core::samplers;
//using namespace transmission_nets::core::computation;
using namespace transmission_nets;

TEST(ContinuousRandomWalkTest, NormalTest) {
    constexpr double TEST_MEAN = 5;
    constexpr double TEST_VARIANCE = 1;
    constexpr int TOTAL_DATA_POINTS = 100;

    struct NormalTestTarget : public core::computation::Computation<core::computation::Likelihood>,
                              public core::abstract::Observable<NormalTestTarget>,
                              public core::abstract::Cacheable<NormalTestTarget>,
                              public core::abstract::Checkpointable<NormalTestTarget, core::computation::Likelihood> {

        explicit NormalTestTarget(core::parameters::Parameter<double> &mean) : mean_(mean) {
            boost::random::normal_distribution<> dist{0, 1};
            for (int i = 0; i < TOTAL_DATA_POINTS; ++i) {
                data_(i) = dist(r) * TEST_VARIANCE + TEST_MEAN;
            }

            mean_.registerCacheableCheckpointTarget(this);
            mean_.add_post_change_listener([=, this]() {
                dirty = true;
            });
            value_ = calculateValue(data_, mean_.value());
        };

        [[nodiscard]] bool isDirty() const {
            return dirty;
        }

        core::computation::Likelihood value() {
            if (dirty) {
                value_ = calculateValue(data_, mean_.value());
                dirty = false;
                std::cout << value_ << " " << mean_.value() << std::endl;
            }

            return value_;
        }

        [[nodiscard]] static core::computation::Likelihood calculateValue(const Eigen::Array<double, TOTAL_DATA_POINTS, 1> &data, double mean) {
            boost::math::normal d(mean, TEST_VARIANCE);
            core::computation::Likelihood val = 0;
            for (auto i = 0; i < data.size(); ++i) {
                val += log(boost::math::pdf(d, data[i]));
            }
            return val;
        }

        boost::random::mt19937 r;
        core::parameters::Parameter<double> &mean_;
        Eigen::Array<double, TOTAL_DATA_POINTS, 1> data_;
        core::computation::Likelihood value_;
        bool dirty{true};
    };


    core::parameters::Parameter<double> myMean(30);
    NormalTestTarget myTestTar(myMean);
    boost::random::mt19937 r;


    core::samplers::ContinuousRandomWalk sampler(myMean, myTestTar, &r, 3);
    sampler.setAdaptationRate(2);

    int i = 20000;
    while (i > 0) {
        i--;
        sampler.update();
        sampler.adapt();
        std::cout << "Variance:" << sampler.variance() << std::endl;
    }

    Eigen::Array<double, 5000, 1> results;
    i = 5000;
    while (i > 0) {
        i--;
        sampler.update();
        results(i) = myMean.value();
    }

    auto resultsMean = results.mean();
    auto resultsStdDev = std::sqrt((results - results.mean()).square().sum() / results.size());
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Variance: " << sampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << sampler.acceptanceRate() << std::endl;
    EXPECT_GE(resultsStdDev * 3, resultsMean - TEST_MEAN);
    EXPECT_LE(-resultsStdDev * 3, resultsMean - TEST_MEAN);
}

TEST(ContinuousRandomWalkTest, DoubleWellTest) {
    struct DoubleWellTestTarget :  public core::computation::Computation<core::computation::Likelihood>,
                                   public core::abstract::Observable<DoubleWellTestTarget>,
                                   public core::abstract::Cacheable<DoubleWellTestTarget>,
                                   public core::abstract::Checkpointable<DoubleWellTestTarget, core::computation::Likelihood> {

        DoubleWellTestTarget(core::parameters::Parameter<double> &x, core::parameters::Parameter<double> &y) : x_(x), y_(y) {
            x_.registerCacheableCheckpointTarget(this);
            x_.add_post_change_listener([=, this]() {
                this->dirty = true;
            });

            y_.registerCacheableCheckpointTarget(this);
            y_.add_post_change_listener([=, this]() {
                this->dirty = true;
            });

            this->value_ = calculateValue(x_.value(), y_.value());
        };

        core::computation::Likelihood value() override {
            if (dirty) {
                value_ = calculateValue(x_.value(), y_.value());
                dirty = false;
            }
            return value_;
        }

        [[nodiscard]] double calculateValue(double x, double y) const {
            return -((.25 * a_ * std::pow(x, 4)) -
                     (.5 * b_ * std::pow(x, 2)) +
                     (c_ * x) +
                     (.5 * d_ * std::pow(y, 2)));
        }

        [[nodiscard]] bool isDirty() const {
            return dirty;
        }

        core::parameters::Parameter<double> &x_;
        core::parameters::Parameter<double> &y_;
        double a_ = 1;
        double b_ = 6;
        double c_ = 1;
        double d_ = 1;
        bool dirty = true;
        core::computation::Likelihood value_;
    };

    core::parameters::Parameter<double> x(0);
    core::parameters::Parameter<double> y(0);
    DoubleWellTestTarget myTestTar(x, y);
    boost::random::mt19937 r;

    core::samplers::ContinuousRandomWalk xSampler(x, myTestTar, &r, 10, .1, 100);
    core::samplers::ContinuousRandomWalk ySampler(y, myTestTar, &r, 10, .1, 100);

    int i = 20000;
    while (i > 0) {
        i--;
        xSampler.update();
        xSampler.adapt();
        ySampler.update();
        ySampler.adapt();
    }

    constexpr int total_samples = 2000;
    Eigen::Array<double, total_samples, 1> xResults;
    Eigen::Array<double, total_samples, 1> yResults;
    i = total_samples;
    while (i > 0) {
        i--;
        xSampler.update();
        ySampler.update();
        xResults(i) = x.value();
        yResults(i) = y.value();
    }

    auto resultsMean = xResults.mean();
    auto resultsStdDev = std::sqrt((xResults - xResults.mean()).square().sum() / xResults.size());
    std::cout << "Mean: " << resultsMean << std::endl;
    std::cout << "StdDev: " << resultsStdDev << std::endl;
    std::cout << "Variance: " << xSampler.variance() << std::endl;
    std::cout << "Acceptance Ratio: " << xSampler.acceptanceRate() << std::endl;
}