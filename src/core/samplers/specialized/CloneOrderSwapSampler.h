//
// Moves an infection to the other side of an infection with the same genome, exchanging their latent parents.
//

#ifndef TRANSMISSION_NETWORKS_APP_CLONEORDERSWAPSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_CLONEORDERSWAPSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"
#include "core/samplers/genetics/GenomePartition.h"
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
     * A Metropolis-Hastings move on an infection's duration and on its and a partner's latent parents. The partner is
     * chosen uniformly among the infection's clones: the infections that can infect it and be infected by it and have
     * its latent genome. The duration is drawn from `durationShape` (the duration prior, which rounds the duration to
     * whole days) restricted to the durations that put the infection on the other side of the partner in the
     * ordering, and the two latent parents' genomes are exchanged. Genomes are unchanged, so the reverse move picks
     * the same partner with the same probability; the correction is the ratio of the duration densities.
     *
     * Of two clones, whichever came first has a latent parent that carries its genome, while the other's latent
     * parent is barely constrained and drifts. Changing which came first one move at a time then needs the drifted
     * latent parent to carry the genome again before the order can change, and in a clonal cluster chains stay with
     * whichever order they reached first.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    class CloneOrderSwapSampler : public AbstractSampler {
    public:
        struct Partner {
            std::shared_ptr<InfectionEventImpl> infection;
            std::shared_ptr<InfectionEventImpl> latentParent;
        };

        CloneOrderSwapSampler(std::shared_ptr<InfectionEventImpl> infection,
                              std::shared_ptr<InfectionEventImpl> latentParent,
                              std::vector<Partner> partners,
                              std::shared_ptr<T> target,
                              std::shared_ptr<Engine> rng,
                              double minDuration,
                              double maxDuration,
                              std::vector<double> durationShape);

        // The infections that can infect infections[i] and be infected by it, with their latent parents (latentParents
        // is parallel to infections).
        template<typename AllowedRelationships>
        static std::vector<Partner> mutualPartners(const std::vector<std::shared_ptr<InfectionEventImpl>>& infections,
                                                   const std::vector<std::shared_ptr<InfectionEventImpl>>& latentParents,
                                                   AllowedRelationships& allowed,
                                                   size_t i);

        void update() noexcept override;

        [[nodiscard]] unsigned int acceptances() const noexcept { return acceptances_; }
        [[nodiscard]] unsigned int rejections() const noexcept { return rejections_; }

    private:
        // The durations [first, second) for which the infection was infected before the partner (`before`) or after it.
        [[nodiscard]] std::pair<double, double> window(const InfectionEventImpl& partner, bool before) const;

        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<InfectionEventImpl> latentParent_;
        std::vector<Partner> partners_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double minDuration_;
        double maxDuration_;
        DurationShape durationShape_;

        boost::random::uniform_01<> uniform_{};
        unsigned int acceptances_ = 0;
        unsigned int rejections_  = 0;

        std::vector<size_t> clones_{};// indices into partners_
        std::vector<std::shared_ptr<typename InfectionEventImpl::LocusType>> exchanged_{};
    };

    template<typename T, typename Engine, typename InfectionEventImpl>
    CloneOrderSwapSampler<T, Engine, InfectionEventImpl>::CloneOrderSwapSampler(std::shared_ptr<InfectionEventImpl> infection,
                                                                               std::shared_ptr<InfectionEventImpl> latentParent,
                                                                               std::vector<Partner> partners,
                                                                               std::shared_ptr<T> target,
                                                                               std::shared_ptr<Engine> rng,
                                                                               const double minDuration,
                                                                               const double maxDuration,
                                                                               std::vector<double> durationShape)
        : infection_(std::move(infection)),
          latentParent_(std::move(latentParent)),
          partners_(std::move(partners)),
          target_(std::move(target)),
          rng_(std::move(rng)),
          minDuration_(minDuration),
          maxDuration_(maxDuration),
          durationShape_(std::move(durationShape)) {}

    template<typename T, typename Engine, typename InfectionEventImpl>
    template<typename AllowedRelationships>
    std::vector<typename CloneOrderSwapSampler<T, Engine, InfectionEventImpl>::Partner> CloneOrderSwapSampler<T, Engine, InfectionEventImpl>::mutualPartners(
            const std::vector<std::shared_ptr<InfectionEventImpl>>& infections,
            const std::vector<std::shared_ptr<InfectionEventImpl>>& latentParents,
            AllowedRelationships& allowed,
            const size_t i) {
        std::vector<Partner> out{};
        const auto parents = allowed.allowedParents(infections[i]);
        for (size_t j = 0; j < infections.size(); ++j) {
            if (std::ranges::find(parents, infections[j]) == parents.end()) {
                continue;
            }
            if (const auto theirs = allowed.allowedParents(infections[j]); std::ranges::find(theirs, infections[i]) != theirs.end()) {
                out.push_back({infections[j], latentParents[j]});
            }
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    std::pair<double, double> CloneOrderSwapSampler<T, Engine, InfectionEventImpl>::window(const InfectionEventImpl& partner, const bool before) const {
        // Infected before the partner means observation time - duration < the partner's infection time.
        const double boundary = infection_->observationTime()->value() - partner.infectionTime();
        return before ? std::pair{std::max(minDuration_, boundary), maxDuration_} : std::pair{minDuration_, std::min(maxDuration_, boundary)};
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void CloneOrderSwapSampler<T, Engine, InfectionEventImpl>::update() noexcept {
        const SAMPLER_STATE_ID stateId = CloneOrderSwapID;
        clones_.clear();
        for (size_t p = 0; p < partners_.size(); ++p) {
            if (genetics::GenomePartition<InfectionEventImpl>::sameGenome(*partners_[p].infection, *infection_)) {
                clones_.push_back(p);
            }
        }
        if (clones_.empty()) {
            return;
        }
        boost::random::uniform_int_distribution<size_t> pickClone(0, clones_.size() - 1);
        const auto& partner = partners_[clones_[pickClone(*rng_)]];

        auto duration              = infection_->infectionDuration();
        const double current       = duration->value();
        const bool before          = infection_->infectionTime() < partner.infection->infectionTime();
        const auto [first, second] = window(*partner.infection, not before);
        const double proposedMass  = second > first ? durationShape_.mass(first, second) : 0.0;
        if (proposedMass <= 0) {
            return;// the infection cannot be on the other side
        }
        // The reverse move draws the current duration from the current side.
        const auto [backFirst, backSecond] = window(*partner.infection, before);
        const double currentDensity        = durationShape_.density(current);
        const double currentLogQ           = currentDensity > 0 ? std::log(currentDensity / durationShape_.mass(backFirst, backSecond)) : -std::numeric_limits<double>::infinity();
        const double currentLik            = target_->value();

        const double proposed = durationShape_.draw(first, second, proposedMass, *rng_);
        duration->saveState(stateId);
        duration->setValue(proposed);
        exchanged_.clear();
        for (const auto& locus : latentParent_->loci()) {
            auto mine   = latentParent_->latentGenotype(locus);
            auto theirs = partner.latentParent->latentGenotype(locus);
            if (mine->value() != theirs->value()) {
                mine->saveState(stateId);
                theirs->saveState(stateId);
                const auto value = mine->value();
                mine->setValue(theirs->value());
                theirs->setValue(value);
                exchanged_.push_back(locus);
            }
        }

        const double proposedLogQ    = std::log(durationShape_.density(proposed) / proposedMass);
        const double acceptanceRatio = target_->value() - currentLik + currentLogQ - proposedLogQ;
        if (std::log(uniform_(*rng_)) <= acceptanceRatio) {
            ++acceptances_;
            countProposal(true);
            duration->acceptState();
            for (const auto& locus : exchanged_) {
                latentParent_->latentGenotype(locus)->acceptState();
                partner.latentParent->latentGenotype(locus)->acceptState();
            }
        } else {
            ++rejections_;
            countProposal(false);
            duration->restoreState(stateId);
            for (const auto& locus : exchanged_) {
                latentParent_->latentGenotype(locus)->restoreState(stateId);
                partner.latentParent->latentGenotype(locus)->restoreState(stateId);
            }
        }
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_CLONEORDERSWAPSAMPLER_H
