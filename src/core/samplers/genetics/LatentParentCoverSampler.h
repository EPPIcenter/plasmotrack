//
// Redraws an infection's latent parent close to the infection.
//

#ifndef TRANSMISSION_NETWORKS_APP_LATENTPARENTCOVERSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_LATENTPARENTCOVERSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"

#include <boost/random.hpp>

#include <cmath>
#include <memory>

namespace transmission_nets::core::samplers::genetics {

    /*
     * A Metropolis-Hastings move on all of a latent parent's genotypes at once, proposed independently of their current
     * values: at each locus every allele the infection carries is present with probability 1 - rate and every other
     * allele with probability rate, redrawn until some allele is present.
     *
     * While an infection has a sampled parent its latent parent is barely constrained and drifts, and the infection
     * can then only lose that parent if its latent parent happens to carry all of its alleles again. This move keeps
     * a latent parent that could have infected it within reach. Every genotype can be proposed, so it can leave a
     * latent parent that has drifted. Genotypes with more than maxCOI alleles at a locus have probability zero and
     * are rejected without evaluating the model.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    class LatentParentCoverSampler : public AbstractSampler {
        using Genetics = typename InfectionEventImpl::GeneticsType;

    public:
        LatentParentCoverSampler(std::shared_ptr<InfectionEventImpl> infection, std::shared_ptr<InfectionEventImpl> latentParent, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int maxCOI, double rate = 0.01)
            : infection_(std::move(infection)), latentParent_(std::move(latentParent)), target_(std::move(target)), rng_(std::move(rng)), maxCOI_(maxCOI), rate_(rate) {}

        void update() noexcept override;

        // Log density of the proposal at the latent parent's current genotypes.
        [[nodiscard]] double logProposalDensity() const;

        [[nodiscard]] unsigned int acceptances() const noexcept { return acceptances_; }
        [[nodiscard]] unsigned int rejections() const noexcept { return rejections_; }

    private:
        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<InfectionEventImpl> latentParent_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int maxCOI_;
        double rate_;

        boost::random::uniform_01<> uniform_{};
        unsigned int acceptances_ = 0;
        unsigned int rejections_  = 0;
    };

    template<typename T, typename Engine, typename InfectionEventImpl>
    double LatentParentCoverSampler<T, Engine, InfectionEventImpl>::logProposalDensity() const {
        double total = 0;
        for (const auto& locus : latentParent_->loci()) {
            const auto& latent = latentParent_->latentGenotype(locus)->value();
            const auto& child  = infection_->latentGenotype(locus)->value();
            const auto kept    = Genetics::truePositiveCount(child, latent); // the child's alleles present
            const auto dropped = Genetics::falseNegativeCount(child, latent);// the child's alleles absent
            const auto added   = Genetics::falsePositiveCount(child, latent);// other alleles present
            const auto absent  = Genetics::trueNegativeCount(child, latent); // other alleles absent
            const double allAbsent = std::pow(rate_, child.totalPositiveCount()) * std::pow(1.0 - rate_, child.totalAlleles() - child.totalPositiveCount());
            total += kept * std::log(1.0 - rate_) + dropped * std::log(rate_) + added * std::log(rate_) + absent * std::log(1.0 - rate_) - std::log1p(-allAbsent);
        }
        return total;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void LatentParentCoverSampler<T, Engine, InfectionEventImpl>::update() noexcept {
        const SAMPLER_STATE_ID stateId = LatentParentCoverID;
        const double currentLik        = target_->value();
        const double currentLogQ       = logProposalDensity();

        bool impossible = false;
        for (const auto& locus : latentParent_->loci()) {
            auto genotype     = latentParent_->latentGenotype(locus);
            const auto& child = infection_->latentGenotype(locus)->value();
            genotype->saveState(stateId);
            auto proposal = genotype->value();
            do {
                proposal.reset();
                for (unsigned int i = 0; i < proposal.totalAlleles(); ++i) {
                    proposal.set(i, uniform_(*rng_) < (child.allele(i) ? 1.0 - rate_ : rate_));
                }
            } while (proposal.totalPositiveCount() == 0);
            impossible = impossible or proposal.totalPositiveCount() > maxCOI_;
            if (proposal != genotype->value()) {
                genotype->setValue(proposal);
            }
        }
        if (impossible) {
            ++rejections_;
            countProposal(false);
            for (const auto& locus : latentParent_->loci()) {
                latentParent_->latentGenotype(locus)->restoreState(stateId);
            }
            return;
        }

        const double proposedLogQ    = logProposalDensity();
        const double acceptanceRatio = target_->value() - currentLik + currentLogQ - proposedLogQ;
        if (std::log(uniform_(*rng_)) <= acceptanceRatio) {
            ++acceptances_;
            countProposal(true);
            for (const auto& locus : latentParent_->loci()) {
                latentParent_->latentGenotype(locus)->acceptState();
            }
        } else {
            ++rejections_;
            countProposal(false);
            for (const auto& locus : latentParent_->loci()) {
                latentParent_->latentGenotype(locus)->restoreState(stateId);
            }
        }
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_LATENTPARENTCOVERSAMPLER_H
