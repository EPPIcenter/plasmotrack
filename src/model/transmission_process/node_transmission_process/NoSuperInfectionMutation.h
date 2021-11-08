//
// Created by Maxwell Murphy on 6/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H

#include <utility>

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include "core/containers/ParentSet.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/datatypes/Matrix.h"

namespace transmission_nets::model::transmission_process {

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    class NoSuperInfectionMutation : public core::computation::Computation<core::datatypes::LogProbabilityTransitionMatrix<2>>,
        public core::abstract::Observable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>,
        public core::abstract::Cacheable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>,
        public core::abstract::Checkpointable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>, core::datatypes::LogProbabilityTransitionMatrix<2>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
    public:
        explicit NoSuperInfectionMutation(p_ParameterDouble mutProb, p_ParameterDouble lossProb, std::shared_ptr<InterTransmissionProbImpl> intp);

        core::datatypes::LogProbabilityTransitionMatrix<2> value() noexcept override;

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        Likelihood peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        Likelihood calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        Likelihood peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);


    private:
        friend class core::abstract::Checkpointable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>, core::datatypes::LogProbabilityTransitionMatrix<2>>;
        friend class core::abstract::Cacheable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>;

        p_ParameterDouble mutProb_;
        p_ParameterDouble lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> intp_;
    };


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::NoSuperInfectionMutation(
            p_ParameterDouble mutProb, p_ParameterDouble lossProb, std::shared_ptr<InterTransmissionProbImpl> intp) : mutProb_(std::move(mutProb)), lossProb_(std::move(lossProb)), intp_(std::move(intp)) {
        mutProb_->registerCacheableCheckpointTarget(this);
        mutProb_->add_post_change_listener([=, this]() { this->setDirty(); });

        lossProb_->registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([=, this]() { this->setDirty(); });

        intp_->registerCacheableCheckpointTarget(this);
        intp_->add_set_dirty_listener([=, this]() { this->setDirty(); });

        value_.setZero();
        this->setDirty();
        this->value();
    }



    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    core::datatypes::LogProbabilityTransitionMatrix<2> NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            core::datatypes::SquareMatrix<double, 2> t_mat;
            t_mat <<    1 - mutProb_->value(),   mutProb_->value(),
                        lossProb_->value(),      1 - lossProb_->value();

            auto tmp = t_mat;
            value_ = tmp * intp_->value()(1);

            for (int i = 2; i <= MaxTransmissions; ++i) {
                tmp = tmp * t_mat;
                value_ += tmp * intp_->value()(i);
            }

            value_ = this->value_.array().log();
            this->setClean();
        }

        return value_;

    }


template<int MaxTransmissions, typename InterTransmissionProbImpl>
template<typename GeneticsImpl>
Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
    if(ps.size() > 1) {
        return -std::numeric_limits<Likelihood>::infinity();
    }
    double llik = 0.0;
    double total_t00 = 0;
    double total_t01 = 0;
    double total_t10 = 0;
    double total_t11 = 0;

    auto const &childGenotypes = child->latentGenotype();
    auto const childGenotypesIter = childGenotypes.begin();
    for (auto const &parent : ps) {
        auto const &parentGenotypes = parent->latentGenotype();
        auto const parentGenotypesIter = parentGenotypes.begin();
        for (size_t i = 0; i < parentGenotypes.size(); ++i) {
            // assume loci are ordered the same
            auto const &parentGenotypeAtLocus = (*(parentGenotypesIter + i)).second;
            auto const &childGenotypeAtLocus = (*(childGenotypesIter + i)).second;
            const unsigned int t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
            const unsigned int t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
            const unsigned int t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
            const unsigned int t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());

            // no mutation
            total_t00 += t00;

            // mutation
            total_t01 += t01;

            // loss
            total_t10 += t10;

            // no loss
            total_t11 += t11;
        }
    }
    llik += total_t00 * value()(0,0);
    llik += total_t01 * value_(0,1); // direct access after ensuring value has been set clean
    llik += total_t10 * value_(1,0);
    llik += total_t11 * value_(1,1);

    return llik;
}

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        if(ps.size() > 1) {
            return -std::numeric_limits<Likelihood>::infinity();
        }
        double llik = 0.0;
        double total_t00 = 0;
        double total_t01 = 0;
        double total_t10 = 0;
        double total_t11 = 0;

        auto const &childGenotypes = child->latentGenotype();
        auto const childGenotypesIter = childGenotypes.begin();
        for (auto const &parent : ps) {
            auto const &parentGenotypes = parent->latentGenotype();

            auto const parentGenotypesIter = parentGenotypes.begin();
            for (size_t i = 0; i < parentGenotypes.size(); ++i) {
                // assume loci are ordered the same
                auto const &parentGenotypeAtLocus = (*(parentGenotypesIter + i)).second;
                auto const &childGenotypeAtLocus = (*(childGenotypesIter + i)).second;
                const unsigned int t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());

                // no mutation
                total_t00 += t00;

                // mutation
                total_t01 += t01;

                // loss
                total_t10 += t10;

                // no loss
                total_t11 += t11;
            }
        }
        llik += total_t00 * peek()(0,0);
        llik += total_t01 * peek()(0,1);
        llik += total_t10 * peek()(1,0);
        llik += total_t11 * peek()(1,1);

        return llik;
    }


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        Likelihood llik = calculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
    }

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        Likelihood llik = peekCalculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
    }

}


#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
