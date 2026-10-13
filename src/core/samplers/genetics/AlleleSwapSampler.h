//
// Swaps a present allele for an absent one in a latent genotype at one locus.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H

#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"

#include <boost/random.hpp>

#include <cmath>
#include <memory>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    /*
     * A Metropolis move on one locus of a latent genotype: a present allele, chosen uniformly, is removed and an absent
     * allele, chosen uniformly, is added. The numbers of present and absent alleles are unchanged, so the reverse move
     * has the same probability and the proposal is symmetric.
     *
     * Moves that add or remove one allele at a time can only change which allele a single-strain locus carries by
     * passing through an empty locus (impossible) or a second allele (a second strain, usually far less probable). At a
     * locus with no informative observation, a wrong allele then stays, and it can rule out the infection's true parent,
     * since a child carries only alleles its parents carry. A swap makes that change in one step.
     */
    template<typename T, typename Engine, typename GeneticsImpl>
    class AlleleSwapSampler : public AbstractSampler {
    public:
        AlleleSwapSampler(std::shared_ptr<parameters::Parameter<GeneticsImpl>> genotype, std::shared_ptr<T> target, std::shared_ptr<Engine> rng)
            : genotype_(std::move(genotype)), target_(std::move(target)), rng_(std::move(rng)) {}

        void update() noexcept override {
            const SAMPLER_STATE_ID stateId = AlleleSwapID;
            const auto& current            = genotype_->value();
            present_.clear();
            absent_.clear();
            for (unsigned int i = 0; i < current.totalAlleles(); ++i) {
                (current.allele(i) ? present_ : absent_).push_back(i);
            }
            if (present_.empty() or absent_.empty()) {
                return;
            }
            boost::random::uniform_int_distribution<size_t> pickPresent(0, present_.size() - 1), pickAbsent(0, absent_.size() - 1);
            auto proposal = current;
            proposal.set(present_[pickPresent(*rng_)], false);
            proposal.set(absent_[pickAbsent(*rng_)], true);

            const double currentLik = target_->value();
            genotype_->saveState(stateId);
            genotype_->setValue(proposal);
            const double proposedLik = target_->value();
            if (std::log(uniform_(*rng_)) <= proposedLik - currentLik) {
                ++acceptances_;
                countProposal(true);
                genotype_->acceptState();
            } else {
                ++rejections_;
                countProposal(false);
                genotype_->restoreState(stateId);
            }
        }

        [[nodiscard]] unsigned int acceptances() const noexcept { return acceptances_; }
        [[nodiscard]] unsigned int rejections() const noexcept { return rejections_; }

    private:
        std::shared_ptr<parameters::Parameter<GeneticsImpl>> genotype_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        boost::random::uniform_01<> uniform_{};
        unsigned int acceptances_ = 0;
        unsigned int rejections_  = 0;
        std::vector<unsigned int> present_{}, absent_{};
    };

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H
