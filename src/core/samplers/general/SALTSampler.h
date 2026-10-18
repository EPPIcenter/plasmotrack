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
#include "core/samplers/general/DirichletIndependenceSampler.h"

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

    /*
     * Optionally, each update also makes a Dirichlet independence move on the whole simplex, fitted to the
     * sampler's own trajectory (enableDirichletMove). The sampler observes the simplex after every update and refits
     * at the end of each fitting window, from that window's observations alone: the mean, and one precision for all
     * frequencies (the median of each frequency's moment estimate). After the last window the fit is fixed. Until the
     * first fit every update runs the SALT sweep; after it, every saltEvery-th update does, and every moveEvery-th
     * update makes the Dirichlet move (the simplex is observed on every update). The proposal is an equal
     * mixture of Dirichlet(alpha) and Dirichlet(widen * alpha), as in DirichletIndependenceSampler.
     *
     * With a block size, a simplex of more frequencies than that is moved in blocks: a uniformly random partition into
     * blocks of blockSize (the last takes the remainder), each proposed in turn by redrawing the block's shares of its
     * current total from the mixture restricted to the block, which is the conditional of the full proposal. An
     * independence proposal's acceptance falls with dimension, so this keeps loci with many alleles moving.
     */
    template<typename T, typename Engine = boost::random::mt19937>
    class SALTSampler : public AbstractSampler {
    public:
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng);
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance);
        SALTSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, double minVariance, double maxVariance);

        [[nodiscard]] double getVariance(int idx) const noexcept;

        void setLowerLimit_(double lim) noexcept;

        // A block size of 0 moves the whole simplex at once.
        void enableDirichletMove(std::vector<int> fitWindows, double widen, int saltEvery, int moveEvery = 1, size_t blockSize = 0);

        // Accepted and proposed Dirichlet moves.
        [[nodiscard]] SamplerCounts dirichletCounts() const noexcept { return dirichletCounts_; }
        // The fitted concentrations; empty before the first fit.
        [[nodiscard]] const std::vector<double>& dirichletAlpha() const noexcept { return alpha_; }

        void update() noexcept override;

        void adapt() noexcept override;

        double acceptanceRate(int idx) noexcept;

        // One proposal scale per Allele frequency.
        [[nodiscard]] AdaptationState adaptationState() const override;

        void restoreAdaptationState(const AdaptationState& state) override;

    private:
        void saltSweep() noexcept;
        void observe() noexcept;
        void dirichletMove() noexcept;
        void wholeDirichletMove() noexcept;
        void dirichletBlockMove(const std::vector<int>& block) noexcept;
        [[nodiscard]] size_t fitSize() const noexcept { return 5 + 3 * variances_.size(); }

        double sampleProposal(double logitCurr, double variance) noexcept;
        double logMetropolisHastingsAdjustment(LogPQ currLogPQ, LogPQ propLogPQ, int k) noexcept;

        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        boost::random::normal_distribution<> normalDist_{0, 1};
        std::vector<double> variances_{};
        std::vector<double> acceptances_{};
        std::vector<double> rejections_{};

        double minVariance_ = .01;
        double maxVariance_ = 1;

        double adaptationRate_       = 1;
        double targetAcceptanceRate_ = .23;

        double lowerLimit_ = .01;

        int totalUpdates_ = 0;

        bool dirichletMove_ = false;
        std::vector<int> fitWindows_{};
        double widen_  = 1;
        int saltEvery_ = 1;
        int moveEvery_ = 1;
        size_t blockSize_ = 0;
        size_t window_ = 0;// the fitting window being observed; fitWindows_.size() once the fit is fixed
        int observed_  = 0;// observations in this window
        std::uint64_t moveUpdates_ = 0;
        std::vector<double> mean_{}, m2_{}, alpha_{}, wide_{}, draw_{}, blockNarrow_{}, blockWide_{}, blockCurrent_{};
        SamplerCounts dirichletCounts_{};
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
    void SALTSampler<T, Engine>::enableDirichletMove(std::vector<int> fitWindows, const double widen, const int saltEvery, const int moveEvery, const size_t blockSize) {
        dirichletMove_ = true;
        fitWindows_    = std::move(fitWindows);
        widen_         = widen;
        saltEvery_     = std::max(1, saltEvery);
        moveEvery_     = std::max(1, moveEvery);
        blockSize_     = blockSize;
        mean_.assign(variances_.size(), 0.0);
        m2_.assign(variances_.size(), 0.0);
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::update() noexcept {
        if (not dirichletMove_) {
            saltSweep();
            return;
        }
        if (alpha_.empty() or moveUpdates_ % saltEvery_ == 0) {
            saltSweep();
        }
        observe();
        if (not alpha_.empty() and moveUpdates_ % moveEvery_ == 0) {
            dirichletMove();
        }
        ++moveUpdates_;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::observe() noexcept {
        if (window_ >= fitWindows_.size()) {
            return;
        }
        const auto f = parameter_->value().frequencies();
        ++observed_;
        for (size_t a = 0; a < f.size(); ++a) {
            const double d = f[a] - mean_[a];
            mean_[a] += d / observed_;
            m2_[a] += d * (f[a] - mean_[a]);
        }
        if (observed_ < fitWindows_[window_]) {
            return;
        }
        std::vector<double> precisions;
        for (size_t a = 0; a < f.size(); ++a) {
            const double variance = observed_ > 1 ? m2_[a] / (observed_ - 1) : 0.0;
            if (variance > 0) {
                precisions.push_back(mean_[a] * (1 - mean_[a]) / variance - 1);
            }
        }
        if (not precisions.empty()) {
            std::nth_element(precisions.begin(), precisions.begin() + precisions.size() / 2, precisions.end());
            const double precision = std::clamp(precisions[precisions.size() / 2], 1.0, 1e6);
            alpha_.resize(f.size());
            wide_.resize(f.size());
            for (size_t a = 0; a < f.size(); ++a) {
                alpha_[a] = std::max(mean_[a] * precision, 1e-3);
                wide_[a]  = alpha_[a] * widen_;
            }
        }
        ++window_;
        observed_ = 0;
        std::fill(mean_.begin(), mean_.end(), 0.0);
        std::fill(m2_.begin(), m2_.end(), 0.0);
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::dirichletMove() noexcept {
        const size_t k = alpha_.size();
        if (blockSize_ == 0 or k <= blockSize_) {
            wholeDirichletMove();
            return;
        }
        const auto order    = utils::generators::randomSequence(0, static_cast<int>(k), rng_);
        const size_t blocks = k / blockSize_;
        for (size_t b = 0; b < blocks; ++b) {
            const auto end = b + 1 == blocks ? order.end() : order.begin() + static_cast<long>((b + 1) * blockSize_);
            dirichletBlockMove(std::vector<int>(order.begin() + static_cast<long>(b * blockSize_), end));
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::dirichletBlockMove(const std::vector<int>& block) noexcept {
        auto proposal = this->propose(*target_, *rng_);
        ++dirichletCounts_.proposed;
        const auto current = parameter_->value().frequencies();
        blockNarrow_.clear();
        blockWide_.clear();
        blockCurrent_.clear();
        double total = 0;
        for (const int a : block) {
            blockNarrow_.push_back(alpha_[a]);
            blockWide_.push_back(wide_[a]);
            total += current[a];
        }
        for (const int a : block) {
            blockCurrent_.push_back(current[a] / total);
        }
        if (not dirichlet::drawMixture(blockNarrow_, blockWide_, lowerLimit_ / total, *rng_, draw_)) {
            proposal.reject();
            return;
        }
        auto next = current;
        for (size_t i = 0; i < block.size(); ++i) {
            next[block[i]] = total * draw_[i];
        }
        // The block's total is the same before and after, so the map from shares to frequencies has the same
        // Jacobian both ways and the Hastings ratio is that of the shares.
        proposal.change(*parameter_, datatypes::Simplex::exactly(next));
        if (proposal.metropolis(dirichlet::logMixtureDensity(blockNarrow_, blockWide_, blockCurrent_) - dirichlet::logMixtureDensity(blockNarrow_, blockWide_, draw_))) {
            ++dirichletCounts_.accepted;
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::wholeDirichletMove() noexcept {
        auto proposal = this->propose(*target_, *rng_);
        ++dirichletCounts_.proposed;
        if (not dirichlet::drawMixture(alpha_, wide_, lowerLimit_, *rng_, draw_)) {
            proposal.reject();
            return;
        }
        const auto current = parameter_->value().frequencies();
        proposal.change(*parameter_, datatypes::Simplex::exactly(draw_));
        if (proposal.metropolis(dirichlet::logMixtureDensity(alpha_, wide_, current) - dirichlet::logMixtureDensity(alpha_, wide_, draw_))) {
            ++dirichletCounts_.accepted;
        }
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::saltSweep() noexcept {
        auto indices = utils::generators::randomSequence(0, parameter_->value().totalElements(), rng_);
        for (const auto idx : indices) {
            auto proposal = this->propose(*target_, *rng_);
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
            bool outOfBounds = false;
            for (const auto el : propVec) {
                if (el < lowerLimit_) {
                    outOfBounds = true;
                    break;
                }
            }
            // A rejection, and the update goes on to the next frequency: whether the rest are proposed must not
            // depend on this proposal.
            if (outOfBounds) {
                rejections_.at(idx)++;
                proposal.reject();
                continue;
            }
            currentVal.set(propVec);
            proposal.change(*parameter_, currentVal);
            const double adjRatio = logMetropolisHastingsAdjustment(currLogPQ, propLogPQ, currentVal.totalElements());
            if (proposal.metropolis(adjRatio)) {
                acceptances_.at(idx)++;
            } else {
                rejections_.at(idx)++;
            }
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
        // The proposal moves one frequency p in logit space by symmetric noise and rescales the other k - 1 so they
        // sum to 1 - p'. The map's Jacobian, p' (1 - p')^(k - 1) / (p (1 - p)^(k - 1)), is the correction.
        return (propLogPQ.logP[0] - currLogPQ.logP[0]) + (k - 1) * (propLogPQ.logQ[0] - currLogPQ.logQ[0]);
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
    AdaptationState SALTSampler<T, Engine>::adaptationState() const {
        AdaptationState state{.updates = static_cast<std::uint64_t>(totalUpdates_)};
        for (size_t j = 0; j < variances_.size(); ++j) {
            const auto accepted = static_cast<std::uint64_t>(acceptances_[j]);
            state.scales.push_back({.variance = variances_[j], .counts = {.accepted = accepted, .proposed = accepted + static_cast<std::uint64_t>(rejections_[j])}});
        }
        if (dirichletMove_) {
            // window, observations, move updates, Dirichlet accepted and proposed, then per frequency the window's
            // running mean and sum of squares and the fitted concentration (0 before the first fit).
            state.fit = {static_cast<double>(window_), static_cast<double>(observed_), static_cast<double>(moveUpdates_), static_cast<double>(dirichletCounts_.accepted),
                         static_cast<double>(dirichletCounts_.proposed)};
            state.fit.insert(state.fit.end(), mean_.begin(), mean_.end());
            state.fit.insert(state.fit.end(), m2_.begin(), m2_.end());
            for (size_t j = 0; j < variances_.size(); ++j) {
                state.fit.push_back(alpha_.empty() ? 0.0 : alpha_[j]);
            }
        }
        return state;
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::restoreAdaptationState(const AdaptationState& state) {
        requireScales(state, variances_.size());
        for (size_t j = 0; j < variances_.size(); ++j) {
            const auto& scale = state.scales[j];
            variances_[j]     = scale.variance;
            acceptances_[j]   = static_cast<double>(scale.counts.accepted);
            rejections_[j]    = static_cast<double>(scale.counts.proposed - scale.counts.accepted);
        }
        totalUpdates_ = static_cast<int>(state.updates);
        // A state without a fit (from before the move, or with it off) leaves the move to fit afresh.
        if (not dirichletMove_ or state.fit.empty()) {
            return;
        }
        if (state.fit.size() != fitSize()) {
            throw std::invalid_argument("Dirichlet fit has " + std::to_string(state.fit.size()) + " values, sampler has " + std::to_string(fitSize()));
        }
        const size_t k   = variances_.size();
        window_          = static_cast<size_t>(state.fit[0]);
        observed_        = static_cast<int>(state.fit[1]);
        moveUpdates_     = static_cast<std::uint64_t>(state.fit[2]);
        dirichletCounts_ = {.accepted = static_cast<std::uint64_t>(state.fit[3]), .proposed = static_cast<std::uint64_t>(state.fit[4])};
        mean_.assign(state.fit.begin() + 5, state.fit.begin() + 5 + k);
        m2_.assign(state.fit.begin() + 5 + k, state.fit.begin() + 5 + 2 * k);
        alpha_.assign(state.fit.begin() + 5 + 2 * k, state.fit.end());
        if (alpha_[0] == 0.0) {
            alpha_.clear();
            wide_.clear();
        } else {
            wide_ = alpha_;
            for (auto& a : wide_) {
                a *= widen_;
            }
        }
    }

    template<typename T, typename Engine>
    double SALTSampler<T, Engine>::getVariance(int idx) const noexcept {
        return variances_.at(idx);
    }

    template<typename T, typename Engine>
    void SALTSampler<T, Engine>::setLowerLimit_(double lim) noexcept {
        lowerLimit_ = lim;
    }


}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H
