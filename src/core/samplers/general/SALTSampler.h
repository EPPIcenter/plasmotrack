//
// Created by Maxwell Murphy on 4/7/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H

#include <algorithm>

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>
#include <utility>

#include "core/datatypes/Simplex.h"
#include "core/io/serialize.h"
#include "core/samplers/AbstractSampler.h"

#include "core/utils/generators/RandomSequence.h"
#include "core/utils/numerics.h"


namespace transmission_nets::core::samplers {
    //  @article{
    //      doi:10.1080/00949655.2017.1376063,
    //      author = {Hannah M. Director and James Gattiker and Earl Lawrence and Scott Vander Wiel},
    //      title = {Efficient sampling on the simplex with a self-adjusting logit transform proposal},
    //      journal = {Journal of Statistical Computation and Simulation},
    //      volume = {87},
    //      number = {18},
    //      pages = {3521-3536},
    //      year  = {2017},
    //      publisher = {Taylor & Francis},
    //      doi = {10.1080/00949655.2017.1376063},
    //  }

    using LogPQ = utils::LogPQ;

    template<typename T, typename Engine = boost::random::mt19937>
    class SALTSampler : public AbstractSampler {
    public:
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng);
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance);
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, double minVariance, double maxVariance);

        void setAdaptationRate(double adaptationRate);

        void setTargetAcceptanceRate(double targetAcceptanceRate);

        void setMinVariance(double minVariance);

        void setMaxVariance(double maxVariance);

        [[nodiscard]] double getVariance(int idx) const noexcept;

        void setVariance(int idx, double var) noexcept;

        void setLowerLimit_(double lim) noexcept;

        void update() noexcept override;

        void adapt() noexcept override;

        void adapt(unsigned int idx) noexcept override;

        double acceptanceRate(int idx) noexcept;


    private:
        double sampleProposal(double logitCurr, double variance) noexcept;
        double logMetropolisHastingsAdjustment(LogPQ currLogPQ, LogPQ propLogPQ, int k) noexcept;

        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        boost::random::normal_distribution<> normalDist_{0, 1};
        boost::random::uniform_01<> uniformDist_{};

        std::vector<double> variances_{};
        std::vector<double> acceptances_{};
        std::vector<double> rejections_{};

        double minVariance_ = .01;
        double maxVariance_ = 1;

        double adaptationRate_       = 1;
        double targetAcceptanceRate_ = .23;

        double lowerLimit_ = .0001;

        int totalUpdates_ = 0;
    };

    template<typename T, typename Engine>
    SALTSampler<T, Engine>::SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) : parameter_(std::move(parameter)), target_(target), rng_(rng) {
        for (size_t j = 0; j < parameter_->value().totalElements(); ++j) {
            variances_.push_back(1);
            acceptances_.push_back(0);
            rejections_.push_back(0);
        }
    }

    template<typename T, typename Engine>
    SALTSampler<T, Engine>::SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance) : parameter_(std::move(parameter)), target_(target), rng_(rng) {
        for (size_t j = 0; j < parameter_->value().totalElements(); ++j) {
            variances_.push_back(variance);
            acceptances_.push_back(0);
            rejections_.push_back(0);
        }
    }

    template<typename T, typename Engine>
    SALTSampler<T, Engine>::SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, double minVariance, double maxVariance) : parameter_(std::move(parameter)), target_(target), rng_(rng), minVariance_(minVariance), maxVariance_(maxVariance) {
        for (size_t j = 0; j < parameter_->value().totalElements(); ++j) {
            variances_.push_back(variance);
            acceptances_.push_back(0);
            rejections_.push_back(0);
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::update() noexcept {
        const std::string stateId = "SALT";
        auto indices              = utils::generators::randomSequence(0, parameter_->value().totalElements(), rng_);
        for (const auto idx : indices) {
            Likelihood curLik = target_->value();

            datatypes::Simplex currentVal(parameter_->value());
            std::vector<double> logitCurrVec = utils::logit(currentVal.frequencies());
            double logitCurrVal              = logitCurrVec[idx];
            double logitPropVal              = sampleProposal(logitCurrVal, variances_[idx]);

            auto currLogPQ = LogPQ({logitCurrVal});
            auto propLogPQ = LogPQ({logitPropVal});

            std::vector<double> logitPropVec(logitCurrVec);
            logitPropVec.erase(logitPropVec.begin() + idx);
            double ls    = propLogPQ.logQ[0] - utils::logitSum(logitPropVec);
            logitPropVec = utils::logitScale(logitPropVec, ls);
            logitPropVec.insert(logitPropVec.begin() + idx, logitPropVal);
            auto propVec = utils::expit(logitPropVec);

            // check to make sure proposal is within lower limit bounds
            for (const auto el : propVec) {
                if (el < lowerLimit_) {
                    rejections_.at(idx)++;
                    totalUpdates_++;
                    return;
                }
            }

            currentVal.set(propVec);

            parameter_->saveState(stateId);

            assert(!target_->isDirty());
            parameter_->setValue(currentVal);
            assert(target_->isDirty());
            const Likelihood newLik = target_->value();

            const double adjRatio = logMetropolisHastingsAdjustment(currLogPQ, propLogPQ, currentVal.totalElements());

            const Likelihood acceptanceRatio = newLik - curLik + adjRatio;

            const bool accept = log(uniformDist_(*rng_)) <= acceptanceRatio;

            if (accept) {
                acceptances_.at(idx)++;
                parameter_->acceptState();
            } else {
                rejections_.at(idx)++;
                parameter_->restoreState(stateId);
                assert(!(target_->isDirty()));
                assert(curLik == target_->value());
            }
            assert(!target_->isDirty());
        }

        totalUpdates_++;
    }

    template<typename T, typename Engine>
    double SALTSampler<T, Engine>::sampleProposal(double logitCurr, double variance) noexcept {
        double eps = normalDist_(*rng_) * variance;
        return logitCurr + eps;
    }

    template<typename T, typename Engine>
    double SALTSampler<T, Engine>::logMetropolisHastingsAdjustment(utils::LogPQ currLogPQ, utils::LogPQ propLogPQ, int k) noexcept {
        return (currLogPQ.logP[0] - propLogPQ.logP[0]) + (k - 1) * (currLogPQ.logQ[0] - propLogPQ.logQ[0]);
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::adapt(unsigned int idx) noexcept {
        for (unsigned int j = 0; j < parameter_->value().totalElements(); ++j) {
            variances_.at(j) += (acceptanceRate(j) - targetAcceptanceRate_) / std::pow(idx, adaptationRate_);
            variances_.at(j) = std::clamp(variances_.at(j), minVariance_, maxVariance_);
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::adapt() noexcept {
        for (unsigned int j = 0; j < parameter_->value().totalElements(); ++j) {
            variances_.at(j) += (acceptanceRate(j) - targetAcceptanceRate_) / std::pow(totalUpdates_ + 1, adaptationRate_);
            variances_.at(j) = std::clamp(variances_.at(j), minVariance_, maxVariance_);
        }
    }

    template<typename T, typename Engine>
    double SALTSampler<T, Engine>::acceptanceRate(int idx) noexcept {
        if (acceptances_.at(idx) + rejections_.at(idx) == 0) {
            return 0;
        } else {
            return acceptances_.at(idx) / (acceptances_.at(idx) + rejections_.at(idx));
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setAdaptationRate(double adaptationRate) {
        adaptationRate_ = adaptationRate;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setTargetAcceptanceRate(double targetAcceptanceRate) {
        targetAcceptanceRate_ = targetAcceptanceRate;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setMinVariance(double minVariance) {
        minVariance_ = minVariance;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setMaxVariance(double maxVariance) {
        maxVariance_ = maxVariance;
    }

    template<typename T, typename Engine>
    double SALTSampler<T, Engine>::getVariance(int idx) const noexcept {
        return variances_.at(idx);
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setVariance(int idx, double var) noexcept {
        variances_.at(idx) = var;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setLowerLimit_(double lim) noexcept {
        lowerLimit_ = lim;
    }


}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H
