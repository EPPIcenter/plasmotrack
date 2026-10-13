//
// Moves an infection to a different source in one step.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTRELINKSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_PARENTRELINKSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"
#include "core/samplers/genetics/GenomePartition.h"
#include "core/samplers/genetics/ParentSetGenotypeProposal.h"
#include "core/samplers/specialized/DurationShape.h"
#include "core/utils/numerics.h"

#include <boost/random.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    /*
     * A Metropolis-Hastings move on an infection's duration and latent genotypes together. One of its sources is chosen
     * uniformly among those that can have infected it first: its allowed parents infected before its latest possible
     * infection, and its latent parent (whose genotypes are a fixed anchor: the latent parent is summed out of the
     * model). The duration is drawn from `durationShape` (the duration prior, which rounds
     * the duration to whole days) restricted to the durations for which that source came first, and every locus is
     * drawn from the source's genotype by proposeLocus(). The proposal depends only on the other infections, not on
     * this infection's duration or genotypes, so the correction is the ratio of the mixture's densities at the current
     * and proposed states.
     *
     * Moving an infection to another source otherwise means rewriting its genome at every locus and often changing its
     * place in the ordering, which single-locus moves and random walks on the duration cannot do in one step.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    class ParentRelinkSampler : public AbstractSampler {
        using Genetics = typename InfectionEventImpl::GeneticsType;

    public:
        ParentRelinkSampler(std::shared_ptr<InfectionEventImpl> infection,
                            std::shared_ptr<InfectionEventImpl> latentParent,
                            std::vector<std::shared_ptr<InfectionEventImpl>> allowedParents,
                            std::shared_ptr<T> target,
                            std::shared_ptr<Engine> rng,
                            double minDuration,
                            double maxDuration,
                            std::vector<double> durationShape,
                            double fnRate       = 0.01,
                            double fpRate       = 0.01,
                            unsigned int tries = 1);

        void update() noexcept override;

        // Log density of the proposal at the infection's current duration and genotypes.
        [[nodiscard]] double logProposalDensity() const;

        [[nodiscard]] unsigned int acceptances() const noexcept { return acceptances_; }
        [[nodiscard]] unsigned int rejections() const noexcept { return rejections_; }

    private:
        [[nodiscard]] size_t numSources() const noexcept { return allowedParents_.size() + 1; }
        // Sources are the allowed parents, then the latent parent.
        [[nodiscard]] const InfectionEventImpl& source(const size_t s) const { return s < allowedParents_.size() ? *allowedParents_[s] : *latentParent_; }
        // The durations [first, second) for which the source was infected before the infection; empty if second <= first.
        [[nodiscard]] std::pair<double, double> window(size_t s) const;
        // The sources the duration can be drawn for (positive mass), and their masses.
        [[nodiscard]] std::vector<std::pair<size_t, double>> feasibleSources() const;
        // Log probability of drawing the infection's current genotypes from the source. Valid within a
        // logProposalDensity() call, which gathers what it needs.
        [[nodiscard]] double logGenomeProb(size_t s) const;
        [[nodiscard]] double logGenomeProbFrom(const InfectionEventImpl& source) const;

        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<InfectionEventImpl> latentParent_;
        std::vector<std::shared_ptr<InfectionEventImpl>> allowedParents_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double minDuration_;
        double maxDuration_;
        DurationShape durationShape_;
        double fnRate_;
        double fpRate_;
        unsigned int tries_;

        // One draw of the proposal: a duration and the infection's genotypes by locus, in infection_->loci() order.
        struct Try {
            double duration = 0;
            std::vector<Genetics> genotypes{};
        };
        void drawTry(Try& t);
        // Sets the infection to `t`, evaluates the target and the proposal density there, and restores it (evaluated).
        std::pair<double, double> evaluateTry(const Try& t);
        void multipleTryUpdate();
        std::vector<Try> forward_{}, reference_{};

        boost::random::uniform_01<> uniform_{};
        unsigned int acceptances_ = 0;
        unsigned int rejections_  = 0;

        // The genome term depends on a source only through its genome, and in a clonal cluster most allowed parents
        // share one, so it is evaluated once per class of allowed parents with the same genome.
        mutable genetics::GenomePartition<InfectionEventImpl> sourceClasses_{};
        mutable std::vector<double> classLogProb_{};// NaN until evaluated in the current logProposalDensity() call
        mutable genetics::ProposalTarget<Genetics> genotypes_{};
        mutable std::vector<Genetics> sourceGenotypes_{};
    };

    template<typename T, typename Engine, typename InfectionEventImpl>
    ParentRelinkSampler<T, Engine, InfectionEventImpl>::ParentRelinkSampler(std::shared_ptr<InfectionEventImpl> infection,
                                                                           std::shared_ptr<InfectionEventImpl> latentParent,
                                                                           std::vector<std::shared_ptr<InfectionEventImpl>> allowedParents,
                                                                           std::shared_ptr<T> target,
                                                                           std::shared_ptr<Engine> rng,
                                                                           const double minDuration,
                                                                           const double maxDuration,
                                                                           std::vector<double> durationShape,
                                                                           const double fnRate,
                                                                           const double fpRate,
                                                                           const unsigned int tries)
        : infection_(std::move(infection)),
          latentParent_(std::move(latentParent)),
          allowedParents_(std::move(allowedParents)),
          target_(std::move(target)),
          rng_(std::move(rng)),
          minDuration_(minDuration),
          maxDuration_(maxDuration),
          durationShape_(std::move(durationShape)),
          fnRate_(fnRate),
          fpRate_(fpRate),
          tries_(std::max(1u, tries)) {}

    template<typename T, typename Engine, typename InfectionEventImpl>
    std::pair<double, double> ParentRelinkSampler<T, Engine, InfectionEventImpl>::window(const size_t s) const {
        if (s >= allowedParents_.size()) {
            return {minDuration_, maxDuration_};
        }
        // Infected before the source means observation time - duration > the source's infection time.
        const double latest = infection_->observationTime()->value() - allowedParents_[s]->infectionTime();
        return {minDuration_, std::min(maxDuration_, latest)};
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    std::vector<std::pair<size_t, double>> ParentRelinkSampler<T, Engine, InfectionEventImpl>::feasibleSources() const {
        std::vector<std::pair<size_t, double>> out{};
        for (size_t s = 0; s < numSources(); ++s) {
            const auto [first, second] = window(s);
            if (second > first) {
                if (const double mass = durationShape_.mass(first, second); mass > 0) {
                    out.emplace_back(s, mass);
                }
            }
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double ParentRelinkSampler<T, Engine, InfectionEventImpl>::logGenomeProb(const size_t s) const {
        if (s >= allowedParents_.size()) {
            return logGenomeProbFrom(*latentParent_);
        }
        const auto c = sourceClasses_.classOf(s);
        if (std::isnan(classLogProb_[c])) {
            classLogProb_[c] = logGenomeProbFrom(sourceClasses_.member(c));
        }
        return classLogProb_[c];
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double ParentRelinkSampler<T, Engine, InfectionEventImpl>::logGenomeProbFrom(const InfectionEventImpl& from) const {
        genetics::latentGenotypesAlong(*infection_, from, sourceGenotypes_);
        double total = 0;
        for (size_t l = 0; l < sourceGenotypes_.size(); ++l) {
            const std::array<const Genetics*, 1> parents{&sourceGenotypes_[l]};
            total += genetics::logLocusProposalProb(parents, 1, genotypes_.latent[l], genotypes_.observed[l], genotypes_.presence[l], fpRate_);
            if (total == -std::numeric_limits<double>::infinity()) {
                break;
            }
        }
        return total;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double ParentRelinkSampler<T, Engine, InfectionEventImpl>::logProposalDensity() const {
        const double duration = infection_->infectionDuration()->value();
        const double density  = durationShape_.density(duration);
        if (density <= 0) {
            return -std::numeric_limits<double>::infinity();
        }
        const auto sources = feasibleSources();
        sourceClasses_.update(allowedParents_);
        classLogProb_.assign(sourceClasses_.numClasses(), std::numeric_limits<double>::quiet_NaN());
        genotypes_.gather(*infection_, fnRate_);
        std::vector<double> terms{};
        terms.reserve(sources.size());
        for (const auto& [s, mass] : sources) {
            const auto [first, second] = window(s);
            if (duration < first or duration >= second) {
                continue;
            }
            terms.push_back(std::log(density / mass) + logGenomeProb(s));
        }
        if (terms.empty()) {
            return -std::numeric_limits<double>::infinity();
        }
        const double maxTerm = *std::max_element(terms.begin(), terms.end());
        return utils::logSumExpKnownMax(terms.begin(), terms.end(), maxTerm) - std::log(static_cast<double>(sources.size()));
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void ParentRelinkSampler<T, Engine, InfectionEventImpl>::drawTry(Try& t) {
        const auto sources = feasibleSources();
        boost::random::uniform_int_distribution<size_t> pickSource(0, sources.size() - 1);
        const auto [s, mass]       = sources[pickSource(*rng_)];
        const auto [first, second] = window(s);
        t.duration                 = durationShape_.draw(first, second, mass, *rng_);
        t.genotypes.clear();
        for (const auto& locus : infection_->loci()) {
            const std::array<const Genetics*, 1> parents{&source(s).latentGenotype(locus)->value()};
            if (infection_->observedGenotype().contains(locus)) {
                t.genotypes.push_back(genetics::proposeLocus(parents, 1, infection_->observedGenotype(locus)->value(), *rng_, fnRate_, fpRate_));
            } else {
                auto nothingObserved = infection_->latentGenotype(locus)->value();
                nothingObserved.reset();
                t.genotypes.push_back(genetics::proposeLocus(parents, 1, nothingObserved, *rng_, genetics::UNOBSERVED_LOCUS_PRESENCE, fpRate_));
            }
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    std::pair<double, double> ParentRelinkSampler<T, Engine, InfectionEventImpl>::evaluateTry(const Try& t) {
        const SAMPLER_STATE_ID stateId = ParentRelinkTryID;
        auto duration                  = infection_->infectionDuration();
        duration->saveState(stateId);
        duration->setValue(t.duration);
        size_t l = 0;
        for (const auto& locus : infection_->loci()) {
            auto genotype = infection_->latentGenotype(locus);
            genotype->saveState(stateId);
            if (t.genotypes[l] != genotype->value()) {
                genotype->setValue(t.genotypes[l]);
            }
            ++l;
        }
        const double value = target_->value();
        const double logQ  = logProposalDensity();
        duration->restoreState(stateId);
        for (const auto& locus : infection_->loci()) {
            infection_->latentGenotype(locus)->restoreState(stateId);
        }
        target_->value();
        return {value, logQ};
    }

    /*
     * With more than one try, a multiple-try Metropolis move (Liu, Liang and Wong 2000) with independent proposals:
     * `tries_` draws y_k of the relink proposal q, each weighted by w(y) = pi(y) / q(y); one is chosen in proportion to
     * its weight; tries_ - 1 fresh draws of q and the current state x form the reference set; and the move is accepted
     * with probability min(1, sum_k w(y_k) / (w(x) + sum_k w(x*_k))). This favours the source and duration the target
     * prefers among several, so a relink whose source must be found among several candidates is proposed far more
     * often than by one uniform pick. A current state q cannot produce has infinite weight and is never left this way.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    void ParentRelinkSampler<T, Engine, InfectionEventImpl>::multipleTryUpdate() {
        constexpr double negInf = -std::numeric_limits<double>::infinity();
        auto logSumExp          = [](const std::vector<double>& xs) {
            double m = negInf;
            for (const double x : xs) {
                m = std::max(m, x);
            }
            if (m == negInf or m == std::numeric_limits<double>::infinity()) {
                return m;
            }
            double sum = 0;
            for (const double x : xs) {
                sum += std::exp(x - m);
            }
            return m + std::log(sum);
        };

        const double piX  = target_->value();
        const double logQ = logProposalDensity();

        forward_.resize(tries_);
        std::vector<double> weights(tries_, negInf);
        for (unsigned int k = 0; k < tries_; ++k) {
            drawTry(forward_[k]);
            const auto [value, q] = evaluateTry(forward_[k]);
            weights[k]            = q == negInf ? negInf : value - q;
        }
        const double forward = logSumExp(weights);
        if (forward == negInf) {
            ++rejections_;
            countProposal(false);
            return;
        }
        const double u = uniform_(*rng_);
        size_t chosen  = tries_ - 1;
        double acc     = 0;
        for (unsigned int k = 0; k < tries_; ++k) {
            acc += std::exp(weights[k] - forward);
            if (u < acc) {
                chosen = k;
                break;
            }
        }

        reference_.resize(tries_ - 1);
        std::vector<double> references{logQ == negInf ? std::numeric_limits<double>::infinity() : piX - logQ};
        for (unsigned int k = 0; k + 1 < tries_; ++k) {
            drawTry(reference_[k]);
            const auto [value, q] = evaluateTry(reference_[k]);
            references.push_back(q == negInf ? negInf : value - q);
        }
        const double backward = logSumExp(references);

        if (std::log(uniform_(*rng_)) <= forward - backward) {
            const SAMPLER_STATE_ID stateId = ParentRelinkID;
            auto duration                  = infection_->infectionDuration();
            duration->saveState(stateId);
            duration->setValue(forward_[chosen].duration);
            size_t l = 0;
            for (const auto& locus : infection_->loci()) {
                auto genotype = infection_->latentGenotype(locus);
                genotype->saveState(stateId);
                if (forward_[chosen].genotypes[l] != genotype->value()) {
                    genotype->setValue(forward_[chosen].genotypes[l]);
                }
                ++l;
            }
            target_->value();
            ++acceptances_;
            countProposal(true);
            duration->acceptState();
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->acceptState();
            }
        } else {
            ++rejections_;
            countProposal(false);
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void ParentRelinkSampler<T, Engine, InfectionEventImpl>::update() noexcept {
        if (tries_ > 1) {
            multipleTryUpdate();
            return;
        }
        const SAMPLER_STATE_ID stateId = ParentRelinkID;
        const double currentLik        = target_->value();
        const double currentLogQ       = logProposalDensity();

        const auto sources = feasibleSources();
        boost::random::uniform_int_distribution<size_t> pickSource(0, sources.size() - 1);
        const auto [s, mass]       = sources[pickSource(*rng_)];
        const auto [first, second] = window(s);

        const double proposed = durationShape_.draw(first, second, mass, *rng_);
        auto duration         = infection_->infectionDuration();
        duration->saveState(stateId);
        duration->setValue(proposed);
        for (const auto& locus : infection_->loci()) {
            auto genotype = infection_->latentGenotype(locus);
            genotype->saveState(stateId);
            const std::array<const Genetics*, 1> parents{&source(s).latentGenotype(locus)->value()};
            Genetics proposal;
            if (infection_->observedGenotype().contains(locus)) {
                proposal = genetics::proposeLocus(parents, 1, infection_->observedGenotype(locus)->value(), *rng_, fnRate_, fpRate_);
            } else {
                auto nothingObserved = genotype->value();
                nothingObserved.reset();
                proposal = genetics::proposeLocus(parents, 1, nothingObserved, *rng_, genetics::UNOBSERVED_LOCUS_PRESENCE, fpRate_);
            }
            if (proposal != genotype->value()) {
                genotype->setValue(proposal);
            }
        }

        const double proposedLogQ    = logProposalDensity();
        const double acceptanceRatio = target_->value() - currentLik + currentLogQ - proposedLogQ;
        if (std::log(uniform_(*rng_)) <= acceptanceRatio) {
            ++acceptances_;
            countProposal(true);
            duration->acceptState();
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->acceptState();
            }
        } else {
            ++rejections_;
            countProposal(false);
            duration->restoreState(stateId);
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->restoreState(stateId);
            }
        }
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_PARENTRELINKSAMPLER_H
