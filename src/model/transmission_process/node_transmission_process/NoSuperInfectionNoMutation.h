//
// Created by Maxwell Murphy on 2/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/datatypes/Matrix.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/ParentSet.h"

namespace transmission_nets::model::transmission_process {

    template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
    class NoSuperInfectionNoMutation : public core::computation::Computation<core::datatypes::LogProbabilityTransitionMatrix<MAX_COI + 1>>,
                                       public core::abstract::Observable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                                       public core::abstract::Cacheable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                                       public core::abstract::Checkpointable<
                                               NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>,
                                               core::datatypes::LogProbabilityTransitionMatrix<MAX_COI + 1>> {

    public:
        explicit NoSuperInfectionNoMutation(std::shared_ptr<COITransitionProbImpl> coitp, std::shared_ptr<InterTransmissionProbImpl> intp);

        core::datatypes::LogProbabilityTransitionMatrix<MAX_COI + 1> value() noexcept override;

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        Likelihood peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        Likelihood calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
            Likelihood llik = calculateLogLikelihood(child, ps);
            return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
        }

        template<typename GeneticsImpl>
        Likelihood peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
            Likelihood llik = peekCalculateLogLikelihood(child, ps);
            return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
        }

    private:
        friend class core::abstract::Checkpointable<
                NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>,
                core::datatypes::LogProbabilityTransitionMatrix<MAX_COI + 1>
                >;

        friend class core::abstract::Cacheable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>;

        std::shared_ptr<COITransitionProbImpl> coitp_;
        std::shared_ptr<InterTransmissionProbImpl> intp_;
    };

    template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
    NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::NoSuperInfectionNoMutation(
            std::shared_ptr<COITransitionProbImpl> coitp, std::shared_ptr<InterTransmissionProbImpl> intp)
        : coitp_(std::move(coitp)), intp_(std::move(intp)) {
        coitp_->registerCacheableCheckpointTarget(this);
        coitp_->add_set_dirty_listener([=, this]() { this->setDirty(); });

        intp_->registerCacheableCheckpointTarget(this);
        intp_->add_set_dirty_listener([=, this]() { this->setDirty(); });
        this->setDirty();
        this->value();
    }

    template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
    core::datatypes::LogProbabilityTransitionMatrix<MAX_COI + 1>
    NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            auto tmp = coitp_->value();
            this->value_ = tmp * intp_->value()(1);

            for (int i = 2; i <= MAX_TRANSMISSIONS; ++i) {
                tmp = tmp * coitp_->value();
                this->value_ += tmp * intp_->value()(i);
            }

            this->value_ = this->value_.array().log();
            this->setClean();
        }


        return this->value_;
    }

    template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood
    NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::calculateLogLikelihood(
            std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        assert(ps.size() <= 1);
        Likelihood llik = 0.0;
        auto const &childGenotype = child->latentGenotype();
        for (auto const &parent : ps) {
            auto const &parentGenotypes = parent->latentGenotype();
            for (auto const &[locus, parentGenotypeAtLocus] : parentGenotypes) {
                if (childGenotype.contains(locus)) {
                    auto const &childGenotypeAtLocus = childGenotype.at(locus);
                    const unsigned int parentAlleleCount = parentGenotypeAtLocus->value().totalPositiveCount();
                    const unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(
                            parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
                    llik += value()(parentAlleleCount, retainedAlleleCount);
                }
            }
        }
        return llik;
    }

    template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood
    NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::peekCalculateLogLikelihood(
            std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        assert(ps.size() == 1);
        Likelihood llik = 0.0;
        auto const &childGenotype = child->latentGenotype();
        for (auto const &parent : ps) {
            auto const &parentGenotypes = parent->latentGenotype();
            for (auto const &[locus, parentGenotypeAtLocus] : parentGenotypes) {
                if (childGenotype.contains(locus)) {
                    auto const &childGenotypeAtLocus = childGenotype.at(locus);
                    unsigned int parentAlleleCount = parentGenotypeAtLocus->value().totalPositiveCount();
                    unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus->value(),
                                                                                       childGenotypeAtLocus->value());
                    llik += this->peek()(parentAlleleCount, retainedAlleleCount);
                }
            }
        }

        return llik;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H
