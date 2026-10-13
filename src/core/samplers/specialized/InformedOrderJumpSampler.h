//
// Moves an infection to the other side of one of its partners in the ordering, choosing the partner by the target.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFORMEDORDERJUMPSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_INFORMEDORDERJUMPSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"
#include "core/samplers/specialized/DurationShape.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    /*
     * A multiple-try Metropolis move (Liu, Liang and Wong 2000) on an infection's duration. Each try j belongs to one of
     * the infection's partners (infections that can infect it and be infected by it): its kernel T_j draws the duration
     * from `durationShape` (the duration prior) restricted to the durations that put the infection on the other side of
     * partner j from where it is. One duration is drawn per partner, each is weighted by w_j = pi(y_j) / T_j(x, y_j), and
     * one is chosen in proportion to its weight. From the chosen y, a reference duration is drawn for every other
     * partner by its kernel, the reference for the chosen partner being the current duration x, and the move is
     * accepted with probability min(1, sum_j w_j(y_j) / sum_j w_j(x*_j)), where w_j(x*_j) = pi(x*_j) / T_j(y, x*_j).
     * With the importance weights, this is MTM with lambda_j(x, y) = 1 / (T_j(x, y) T_j(y, x)), which is symmetric, so
     * the move leaves the target invariant. A partner whose other side has no duration-prior mass contributes no try.
     *
     * Moving an infection to the other side of a related infection usually means a large change of duration, which
     * random walks rarely make; choosing the partner by the target finds the side the posterior favours.
     *
     * With more partners than maxTries, each update first draws maxTries of them uniformly without replacement,
     * independently of the state, and makes the move among those; the reverse move draws the same subset with the same
     * probability, so the move still leaves the target invariant. maxTries = 0 uses every partner.
     *
     * Every duration tried is evaluated before it is restored, since restoring a change the target has not evaluated
     * leaves the target's cached value wrong.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    class InformedOrderJumpSampler : public AbstractSampler {
    public:
        InformedOrderJumpSampler(std::shared_ptr<InfectionEventImpl> infection,
                                 std::vector<std::shared_ptr<InfectionEventImpl>> partners,
                                 std::shared_ptr<T> target,
                                 std::shared_ptr<Engine> rng,
                                 double minDuration,
                                 double maxDuration,
                                 std::vector<double> durationShape,
                                 size_t maxTries = 0)
            : infection_(std::move(infection)), partners_(std::move(partners)), target_(std::move(target)), rng_(std::move(rng)),
              minDuration_(minDuration), maxDuration_(maxDuration), durationShape_(std::move(durationShape)), maxTries_(maxTries) {}

        void update() noexcept override;

        [[nodiscard]] unsigned int acceptances() const noexcept { return acceptances_; }
        [[nodiscard]] unsigned int rejections() const noexcept { return rejections_; }

    private:
        static constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

        // The durations [first, second) that put the infection on the other side of partner j from where it is now.
        [[nodiscard]] std::pair<double, double> otherSide(size_t j) const;
        // log target at duration d, leaving the duration and the target as they were (both evaluated).
        double logTargetAt(double d);

        static double logSumExp(const std::vector<double>& xs) {
            double m = NEG_INF;
            for (const double x : xs) {
                m = std::max(m, x);
            }
            if (m == NEG_INF) {
                return m;
            }
            double s = 0;
            for (const double x : xs) {
                s += std::exp(x - m);
            }
            return m + std::log(s);
        }

        std::shared_ptr<InfectionEventImpl> infection_;
        std::vector<std::shared_ptr<InfectionEventImpl>> partners_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double minDuration_;
        double maxDuration_;
        DurationShape durationShape_;
        size_t maxTries_;

        boost::random::uniform_01<> uniform_{};
        unsigned int acceptances_ = 0;
        unsigned int rejections_  = 0;

        std::vector<double> tries_{}, weights_{}, references_{};
        std::vector<size_t> active_{};// the partners tried in this update
    };

    template<typename T, typename Engine, typename InfectionEventImpl>
    std::pair<double, double> InformedOrderJumpSampler<T, Engine, InfectionEventImpl>::otherSide(const size_t j) const {
        // Infected before the partner means observation time - duration < the partner's infection time.
        const double boundary = infection_->observationTime()->value() - partners_[j]->infectionTime();
        const bool before     = infection_->infectionTime() < partners_[j]->infectionTime();
        return before ? std::pair{minDuration_, std::min(maxDuration_, boundary)} : std::pair{std::max(minDuration_, boundary), maxDuration_};
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double InformedOrderJumpSampler<T, Engine, InfectionEventImpl>::logTargetAt(const double d) {
        // Its own state id: the reference tries are evaluated while the move's own change is saved.
        const SAMPLER_STATE_ID stateId = InformedOrderJumpTryID;
        auto duration                  = infection_->infectionDuration();
        duration->saveState(stateId);
        duration->setValue(d);
        const double value = target_->value();
        duration->restoreState(stateId);
        target_->value();
        return value;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void InformedOrderJumpSampler<T, Engine, InfectionEventImpl>::update() noexcept {
        const SAMPLER_STATE_ID stateId = InformedOrderJumpID;
        if (partners_.empty()) {
            return;
        }
        active_.resize(partners_.size());
        for (size_t j = 0; j < active_.size(); ++j) {
            active_[j] = j;
        }
        if (maxTries_ > 0 and active_.size() > maxTries_) {
            // A partial Fisher-Yates shuffle: the first maxTries_ entries are a uniform subset.
            for (size_t j = 0; j < maxTries_; ++j) {
                boost::random::uniform_int_distribution<size_t> pick(j, active_.size() - 1);
                std::swap(active_[j], active_[pick(*rng_)]);
            }
            active_.resize(maxTries_);
        }
        const size_t k = active_.size();
        auto duration    = infection_->infectionDuration();
        const double x   = duration->value();
        const double piX = target_->value();

        // Forward tries.
        tries_.assign(k, 0.0);
        weights_.assign(k, NEG_INF);
        bool any = false;
        for (size_t j = 0; j < k; ++j) {
            const auto [first, second] = otherSide(active_[j]);
            const double mass          = second > first ? durationShape_.mass(first, second) : 0.0;
            if (mass <= 0) {
                continue;
            }
            const double y    = durationShape_.draw(first, second, mass, *rng_);
            const double logT = std::log(durationShape_.density(y) / mass);
            tries_[j]         = y;
            weights_[j]       = logTargetAt(y) - logT;
            any               = true;
        }
        const double forward = logSumExp(weights_);
        if (not any or forward == NEG_INF) {
            if (any) {
                ++rejections_;
                countProposal(false);
            }
            return;
        }
        const double u = uniform_(*rng_);
        size_t chosen = k;
        {
            double acc = 0;
            for (size_t j = 0; j < k; ++j) {
                acc += std::exp(weights_[j] - forward);
                if (u < acc) {
                    chosen = j;
                    break;
                }
            }
            if (chosen == k) {// rounding: the last try with any weight
                for (size_t j = k; j-- > 0;) {
                    if (weights_[j] > NEG_INF) {
                        chosen = j;
                        break;
                    }
                }
            }
        }

        // Move to the chosen duration, then draw the reference tries from there.
        duration->saveState(stateId);
        duration->setValue(tries_[chosen]);
        target_->value();
        references_.assign(k, NEG_INF);
        for (size_t j = 0; j < k; ++j) {
            const auto [first, second] = otherSide(active_[j]);
            const double mass          = second > first ? durationShape_.mass(first, second) : 0.0;
            if (j == chosen) {
                // The reference for the chosen try is the current duration, on the side the move left.
                const double density = durationShape_.density(x);
                references_[j]       = mass > 0 and density > 0 ? piX - std::log(density / mass) : NEG_INF;
                continue;
            }
            if (mass <= 0) {
                continue;
            }
            const double r    = durationShape_.draw(first, second, mass, *rng_);
            const double logT = std::log(durationShape_.density(r) / mass);
            references_[j]    = logTargetAt(r) - logT;
        }
        const double backward = logSumExp(references_);

        if (std::log(uniform_(*rng_)) <= forward - backward) {
            ++acceptances_;
            countProposal(true);
            duration->acceptState();
        } else {
            ++rejections_;
            countProposal(false);
            duration->restoreState(stateId);
            target_->value();
        }
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_INFORMEDORDERJUMPSAMPLER_H
