//
// Moves every frequency of a simplex in one Metropolis-Hastings step, in log-ratio space, with a proposal covariance
// learned from the chain while it adapts.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEXLOGRATIOSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEXLOGRATIOSAMPLER_H

#include "core/computation/PartialLikelihood.h"
#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <boost/random.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

namespace transmission_nets::core::samplers::general {

    /*
     * A random walk on a simplex π in centred log-ratio coordinates y = log π - mean(log π). A proposal adds Gaussian
     * noise δ to log π and renormalises: π' ∝ π e^δ. Only δ's part orthogonal to the all-ones vector moves y, and that
     * part is symmetric, so the Hastings ratio is the Jacobian of y → π, Π π'_i / Π π_i. Proposals with a frequency
     * below the lower limit are rejected before any state changes.
     *
     * δ ~ N(0, s² Σ). Σ starts as the identity; while the sampler adapts it becomes the chain's covariance of y, scaled
     * to average variance 1 (plus a small ridge), once there are enough draws, and s is moved by Robbins-Monro towards
     * the target acceptance rate, starting from 2.38 / √(k - 1) times the initial scale. Neither changes outside adaptation, so the chain samples
     * the target once adaptation ends. One proposal evaluates the target once, however many frequencies there are.
     */
    template<typename T, typename Engine>
    class SimplexLogRatioSampler : public AbstractSampler {
    public:
        SimplexLogRatioSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                               double lowerLimit, unsigned int proposalsPerUpdate = 1, double initialScale = 0.1)
            : parameter_(std::move(parameter)), target_(std::move(target)), rng_(std::move(rng)), lowerLimit_(lowerLimit),
              proposalsPerUpdate_(proposalsPerUpdate) {
            k_ = parameter_->value().totalElements();
            assert(k_ >= 2);
            chol_     = Eigen::MatrixXd::Identity(k_, k_);
            logScale_ = std::log(initialScale * 2.38 / std::sqrt(static_cast<double>(k_ - 1)));
            mean_     = Eigen::VectorXd::Zero(k_);
            m2_       = Eigen::MatrixXd::Zero(k_, k_);
        }

        void update() noexcept override;

        // Called once per adapting iteration: records the current point and moves the scale by the last update's
        // acceptance rate.
        void adapt() noexcept override;
        void adapt([[maybe_unused]] unsigned int idx) noexcept override { adapt(); }

        [[nodiscard]] double acceptanceRate() const noexcept {
            return proposed_ == 0 ? 0.0 : static_cast<double>(accepted_) / static_cast<double>(proposed_);
        }
        [[nodiscard]] double scale() const noexcept { return std::exp(logScale_); }

        // The scale only; the learned covariance is not kept, so a resumed Run starts it again from the identity.
        [[nodiscard]] AdaptationState adaptationState() const override {
            return {.scales = {{.variance = logScale_, .counts = {.accepted = accepted_, .proposed = proposed_}}}, .updates = adaptations_};
        }
        void restoreAdaptationState(const AdaptationState& state) override {
            requireScales(state, 1);
            logScale_    = state.scales[0].variance;
            accepted_    = state.scales[0].counts.accepted;
            proposed_    = state.scales[0].counts.proposed;
            adaptations_ = state.updates;
        }

    private:
        static constexpr double TARGET_ACCEPTANCE = 0.234;
        static constexpr unsigned int REFACTOR_EVERY = 25;

        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double lowerLimit_;
        unsigned int proposalsPerUpdate_;
        unsigned int k_ = 0;

        boost::random::normal_distribution<> normal_{0, 1};
        boost::random::uniform_01<> uniform_{};

        double logScale_;
        Eigen::MatrixXd chol_;// lower Cholesky factor of Σ
        // Welford moments of y over adapting iterations.
        Eigen::VectorXd mean_;
        Eigen::MatrixXd m2_;
        std::uint64_t observations_ = 0;

        std::uint64_t accepted_    = 0;
        std::uint64_t proposed_    = 0;
        std::uint64_t adaptations_ = 0;
        unsigned int lastAccepted_ = 0;
        unsigned int lastProposed_ = 0;
    };

    template<typename T, typename Engine>
    void SimplexLogRatioSampler<T, Engine>::update() noexcept {
        const SAMPLER_STATE_ID stateId = SimplexLogRatioID;
        lastAccepted_                  = 0;
        lastProposed_                  = 0;
        std::vector<double> proposal(k_);
        Eigen::VectorXd noise(k_);
        for (unsigned int t = 0; t < proposalsPerUpdate_; ++t) {
            const auto current = parameter_->value().frequencies();
            for (unsigned int i = 0; i < k_; ++i) {
                noise[i] = normal_(*rng_);
            }
            const Eigen::VectorXd delta = std::exp(logScale_) * (chol_ * noise);
            double maxLog               = -std::numeric_limits<double>::infinity();
            for (unsigned int i = 0; i < k_; ++i) {
                proposal[i] = std::log(current[i]) + delta[i];
                maxLog      = std::max(maxLog, proposal[i]);
            }
            double sum = 0;
            for (unsigned int i = 0; i < k_; ++i) {
                proposal[i] = std::exp(proposal[i] - maxLog);
                sum += proposal[i];
            }
            bool outOfBounds = false;
            double logJacobian = 0;
            for (unsigned int i = 0; i < k_; ++i) {
                proposal[i] /= sum;
                outOfBounds = outOfBounds or proposal[i] < lowerLimit_;
                logJacobian += std::log(proposal[i]) - std::log(current[i]);
            }
            ++proposed_;
            ++lastProposed_;
            if (outOfBounds) {
                countProposal(false);
                continue;
            }

            const double currentLik = target_->value();
            parameter_->saveState(stateId);
            parameter_->setValue(datatypes::Simplex(proposal));
            const double acceptanceRatio = target_->value() - currentLik + logJacobian;
            if (std::log(uniform_(*rng_)) <= acceptanceRatio) {
                ++accepted_;
                ++lastAccepted_;
                countProposal(true);
                parameter_->acceptState();
            } else {
                countProposal(false);
                parameter_->restoreState(stateId);
            }
        }
    }

    template<typename T, typename Engine>
    void SimplexLogRatioSampler<T, Engine>::adapt() noexcept {
        ++adaptations_;
        if (lastProposed_ > 0) {
            const double rate = static_cast<double>(lastAccepted_) / static_cast<double>(lastProposed_);
            logScale_ += (rate - TARGET_ACCEPTANCE) / std::pow(static_cast<double>(adaptations_), 0.6);
        }

        const auto freqs = parameter_->value().frequencies();
        Eigen::VectorXd y(k_);
        double meanLog = 0;
        for (unsigned int i = 0; i < k_; ++i) {
            y[i] = std::log(freqs[i]);
            meanLog += y[i];
        }
        y.array() -= meanLog / k_;
        ++observations_;
        const Eigen::VectorXd d = y - mean_;
        mean_ += d / static_cast<double>(observations_);
        m2_ += d * (y - mean_).transpose();

        // Use the learned covariance's shape once there are a few draws per dimension. It is scaled to average variance
        // 1, as the identity has, so that s keeps carrying the step size its tuning reached.
        if (observations_ >= 2 * k_ and observations_ % REFACTOR_EVERY == 0) {
            Eigen::MatrixXd cov = m2_ / static_cast<double>(observations_ - 1);
            const double averageVariance = cov.trace() / k_;
            if (averageVariance > 0) {
                cov /= averageVariance;
                cov.diagonal().array() += 1e-6;
                const Eigen::LLT<Eigen::MatrixXd> llt(cov);
                if (llt.info() == Eigen::Success) {
                    chol_ = llt.matrixL();
                }
            }
        }
    }

}// namespace transmission_nets::core::samplers::general

#endif//TRANSMISSION_NETWORKS_APP_SIMPLEXLOGRATIOSAMPLER_H
