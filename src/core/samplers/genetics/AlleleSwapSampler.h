//
// Swaps a present allele for an absent one in a latent genotype at one locus.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H

#include "core/datatypes/Genotype.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

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
    template<typename T, typename Engine>
    class AlleleSwapSampler : public AbstractSampler {
    public:
        AlleleSwapSampler(std::shared_ptr<parameters::Parameter<core::datatypes::Genotype>> genotype, std::shared_ptr<T> target, std::shared_ptr<Engine> rng)
            : genotype_(std::move(genotype)), target_(std::move(target)), rng_(std::move(rng)) {}

        void update() noexcept override {
            const auto& current = genotype_->value();
            present_.clear();
            absent_.clear();
            for (unsigned int i = 0; i < current.totalAlleles(); ++i) {
                (current.allele(i) ? present_ : absent_).push_back(i);
            }
            if (present_.empty() or absent_.empty()) {
                return;
            }
            boost::random::uniform_int_distribution<size_t> pickPresent(0, present_.size() - 1), pickAbsent(0, absent_.size() - 1);
            auto proposed = current;
            proposed.set(present_[pickPresent(*rng_)], false);
            proposed.set(absent_[pickAbsent(*rng_)], true);

            auto proposal = this->propose(*target_, *rng_);
            proposal.change(*genotype_, proposed);
            proposal.metropolis();
        }
    private:
        std::shared_ptr<parameters::Parameter<core::datatypes::Genotype>> genotype_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        std::vector<unsigned int> present_{}, absent_{};
    };

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_ALLELESWAPSAMPLER_H
