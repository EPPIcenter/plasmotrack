//
// Created by Maxwell Murphy on 6/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H


#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/ParentSet.h"

#include "core/datatypes/Matrix.h"


#include <utility>

namespace transmission_nets::model::transmission_process {
    /*
     * Under a model of no super infection and with mutation, creates a 2x2 transition matrix representing the probability of an
     * allele transitioning state [0,1] x [0,1]
     */
    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    class NoSuperInfectionMutation : public core::computation::Computation<std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>>>,
                                     public core::abstract::Observable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>,
                                     public core::abstract::Cacheable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>,
                                     public core::abstract::Checkpointable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>, std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>>> {

        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;

    public:
        explicit NoSuperInfectionMutation(p_Parameterfloat mutProb, p_Parameterfloat lossProb, std::shared_ptr<InterTransmissionProbImpl> intp);

        std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>> value() noexcept override;

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        Likelihood peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        Likelihood calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);

        template<typename GeneticsImpl>
        Likelihood peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps);


    private:
        friend class core::abstract::Checkpointable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>, std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>>>;
        friend class core::abstract::Cacheable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>;

        p_Parameterfloat mutProb_;
        p_Parameterfloat lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> intp_;
    };


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::NoSuperInfectionMutation(
            p_Parameterfloat mutProb, p_Parameterfloat lossProb, std::shared_ptr<InterTransmissionProbImpl> intp) : mutProb_(std::move(mutProb)), lossProb_(std::move(lossProb)), intp_(std::move(intp)) {
        mutProb_->registerCacheableCheckpointTarget(this);
        mutProb_->add_post_change_listener([=, this]() { this->setDirty(); });

        lossProb_->registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([=, this]() { this->setDirty(); });

        intp_->registerCacheableCheckpointTarget(this);
        intp_->add_set_dirty_listener([=, this]() { this->setDirty(); });

        value_ = std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>>(MaxTransmissions + 1);
        for (auto v : value_) {
            v.setZero();
        }

        this->setDirty();
        this->value();
    }

    /**
     * Return a vector of transition matrices representing the probability of an allele transitioning state [0,1] x [0,1]. Vector is
     * indexed by total number of transmission events between the child and the parent.
     * @return
     */
    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    std::vector<core::datatypes::LogProbabilityTransitionMatrix<2>> NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            core::datatypes::SquareMatrix<float, 2> t_mat;
            t_mat << 1 - mutProb_->value(), mutProb_->value(),
                    lossProb_->value(), 1 - lossProb_->value();

            value_[1] = t_mat.array().log();

            auto tmp = t_mat;

            for (int i = 2; i <= MaxTransmissions; ++i) {
                tmp = tmp * t_mat;
                value_[i] = tmp.array().log();
            }
            this->setClean();
        }
        return value_;
    }


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        if (ps.size() > 1) {
            return -std::numeric_limits<Likelihood>::infinity();
        }
        auto llik_vec = std::vector<Likelihood>(MaxTransmissions + 1, 0);
        auto const& childGenotypes    = child->latentGenotype();
        auto const childGenotypesIter = childGenotypes.begin();
        for (auto const& parent : ps) {
            auto const& parentGenotypes    = parent->latentGenotype();
            auto const parentGenotypesIter = parentGenotypes.begin();
            for (size_t i = 0; i < parentGenotypes.size(); ++i) {
                // assume loci are ordered the same -- should enforce this at construction
                auto const& parentGenotypeAtLocus = (*(parentGenotypesIter + i)).second;
                auto const& childGenotypeAtLocus  = (*(childGenotypesIter + i)).second;
                const unsigned int total_t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
                const unsigned int total_t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
                const unsigned int total_t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());
                const unsigned int total_t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus->value(), childGenotypeAtLocus->value());

                for (int tt = 1; tt <= MaxTransmissions; ++tt) {
                    // tt is the number of transmissions between the child and the parent.
                    // calculate the probability of the genetics given tt transmissions
                    llik_vec[tt] += total_t00 * value()[tt](0, 0) + total_t01 * value()[tt](0, 1) + total_t10 * value()[tt](1, 0) + total_t11 * value()[tt](1, 1);
                }
            }
        }
        for (int tt = 1; tt <= MaxTransmissions; ++tt) {
            // add the log probability of there being tt transmissions
            llik_vec[tt] += std::log(intp_->value()(tt));
        }
        Likelihood llik = core::utils::logSumExp(llik_vec);
        return llik;
    }

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        if (ps.size() > 1) {
            return -std::numeric_limits<Likelihood>::infinity();
        }

        auto llik_vec = std::vector<Likelihood>(MaxTransmissions + 1, 0);
        auto const& childGenotypes    = child->latentGenotype();
        auto const childGenotypesIter = childGenotypes.begin();
        for (auto const& parent : ps) {
            auto const& parentGenotypes = parent->latentGenotype();

            auto const parentGenotypesIter = parentGenotypes.begin();
            for (size_t i = 0; i < parentGenotypes.size(); ++i) {
                // assume loci are ordered the same
                auto const& parentGenotypeAtLocus = (*(parentGenotypesIter + i)).second;
                auto const& childGenotypeAtLocus  = (*(childGenotypesIter + i)).second;
                const unsigned int total_t00            = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int total_t01            = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int total_t10            = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                const unsigned int total_t11            = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus->peek(), childGenotypeAtLocus->peek());
                for (int tt = 1; tt <= MaxTransmissions; ++tt) {
                    // tt is the number of transmissions between the child and the parent.
                    // calculate the probability of the genetics given tt transmissions
                    llik_vec[tt] += total_t00 * peek()[tt](0, 0) + total_t01 * peek()[tt](0, 1) + total_t10 * peek()[tt](1, 0) + total_t11 * peek()[tt](1, 1);
                }
            }
        }
        for (int tt = 1; tt <= MaxTransmissions; ++tt) {
            // add the log probability of there being tt transmissions
            llik_vec[tt] += std::log(intp_->peek()(tt));
        }
        Likelihood llik = core::utils::logSumExp(llik_vec);
        return llik;
    }


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        Likelihood llik = calculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
    }

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& ps) {
        Likelihood llik = peekCalculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<Likelihood>::infinity() ? exp(llik) : 0;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
