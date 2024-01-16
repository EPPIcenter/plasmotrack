//
// Created by mmurphy on 11/8/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_SUPERINFECTIONNOMUTATION_H
#define TRANSMISSION_NETWORKS_APP_SUPERINFECTIONNOMUTATION_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"

#include <cmath>
#include <memory>

namespace transmission_nets::model::transmission_process {

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    class SuperInfectionNoMutation : public core::computation::Computation<float>,
                                     public core::abstract::Observable<SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>,
                                     public core::abstract::Cacheable<SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>,
                                     public core::abstract::Checkpointable<SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>, core::computation::Computation<float>> {


    public:
        explicit SuperInfectionNoMutation(std::shared_ptr<InterTransmissionProbImpl> intp, std::shared_ptr<core::parameters::Parameter<float>> lossProb);

        // logp(y = 0 | y' = 1)
        float value() noexcept override;

        template<typename GeneticsImpl>
        core::computation::Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        core::computation::Likelihood calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        core::computation::Likelihood peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        core::computation::Likelihood peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

    private:
        friend class core::abstract::Checkpointable<SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>, core::computation::Computation<float>>;
        friend class core::abstract::Cacheable<SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>;

        std::shared_ptr<core::parameters::Parameter<float>> lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> intp_;
    };

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::SuperInfectionNoMutation(std::shared_ptr<InterTransmissionProbImpl> intp, std::shared_ptr<core::parameters::Parameter<float>> lossProb) : intp_(std::move(intp)), lossProb_(std::move(lossProb)) {
        lossProb_->registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([=, this]() { this->setDirty(); });

        intp_->registerCacheableCheckpointTarget(this);
        intp_->add_set_dirty_listener([=, this]() { this->setDirty(); });
        this->setDirty();
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    float SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            this->value_ = 0;
            for (int ii = 1; ii <= MAX_TRANSMISSIONS; ++ii) {
                this->value_ += (1 - std::pow(1 - lossProb_->value(), ii)) * intp_->value()(ii);
            }
            this->value_ = std::log(this->value_);
            this->setClean();
        }
        return this->value_;
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    core::computation::Likelihood SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::calculateLogLikelihood(
            std::shared_ptr<core::containers::Infection<GeneticsImpl>> child,
            const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        float llik = 0.0;

        // logp(y = 0 | y' = 1)
        float logProbLost            = this->value();

        auto const& childGenotypes    = child->latentGenotype();
        auto const childGenotypesIter = childGenotypes.begin();
        for (size_t ii = 0; ii < childGenotypes.size(); ++ii) {
            const auto& childGenotype = *(childGenotypesIter + ii).second->value();
            std::vector<int> presenceVector(childGenotype.totalAlleles(), 0);
            // For each allele, check how many parents in the parent set have the allele
            for (auto const& parent : ps) {
                const auto& parentGenotype = *(parent->latentGenotype().begin() + ii).second->value();
                for (int jj = 0; jj < childGenotype.totalAlleles(); ++jj) {
                    presenceVector[jj] += parentGenotype.allele(jj);
                }
            }
            for (int jj = 0; jj < childGenotype.totalAlleles(); ++jj) {
                llik += presenceVector[jj] * (((1 - childGenotype.allele(jj)) * logProbLost) + (childGenotype.allele(jj) * (1 - logProbLost)));
            }
        }

        return llik;
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    core::computation::Likelihood SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::peekCalculateLogLikelihood(
            std::shared_ptr<core::containers::Infection<GeneticsImpl>> child,
            const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        float llik = 0.0;
        // logp(y = 0 | y' = 1)
        float logProbLost            = this->peek();
        auto const& childGenotypes    = child->latentGenotype();
        auto const childGenotypesIter = childGenotypes.begin();
        for (size_t ii = 0; ii < childGenotypes.size(); ++ii) {
            const auto& childGenotype = *(childGenotypesIter + ii).second->peek();
            std::vector<int> presenceVector(childGenotype.totalAlleles(), 0);
            // For each allele, check how many parents in the parent set have the allele
            for (auto const& parent : ps) {
                const auto& parentGenotype = *(parent->latentGenotype().begin() + ii).second->peek();
                for (int jj = 0; jj < childGenotype.totalAlleles(); ++jj) {
                    presenceVector[jj] += parentGenotype.allele(jj);
                }
            }

            for (int jj = 0; jj < childGenotype.totalAlleles(); ++jj) {
                llik += presenceVector[jj] * (((1 - childGenotype.allele(jj)) * logProbLost) + (childGenotype.allele(jj) * (1 - logProbLost)));
            }
        }

        return llik;
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    core::computation::Likelihood SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        core::computation::Likelihood llik = calculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<core::computation::Likelihood>::infinity() ? std::exp(llik) : 0;
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    core::computation::Likelihood SuperInfectionNoMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        core::computation::Likelihood llik = peekCalculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<core::computation::Likelihood>::infinity() ? std::exp(llik) : 0;
    }


}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_SUPERINFECTIONNOMUTATION_H
