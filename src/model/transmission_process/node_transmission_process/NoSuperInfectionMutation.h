//
// Created by Maxwell Murphy on 6/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H

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

    public:
        explicit NoSuperInfectionMutation(core::parameters::Parameter<double> &mutProb, core::parameters::Parameter<double> &lossProb, InterTransmissionProbImpl &intp);

        core::datatypes::LogProbabilityTransitionMatrix<2> value() noexcept override;

        template<typename GeneticsImpl>
        double calculateLogLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        double peekCalculateLogLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        double calculateLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);

        template<typename GeneticsImpl>
        double peekCalculateLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps);


    private:
        friend class core::abstract::Checkpointable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>, core::datatypes::LogProbabilityTransitionMatrix<2>>;
        friend class core::abstract::Cacheable<NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>>;

        core::parameters::Parameter<double> &mutProb_;
        core::parameters::Parameter<double> &lossProb_;
        InterTransmissionProbImpl &intp_;
    };


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::NoSuperInfectionMutation(
            core::parameters::Parameter<double> &mutProb, core::parameters::Parameter<double> &lossProb, InterTransmissionProbImpl &intp):mutProb_(mutProb), lossProb_(lossProb), intp_(intp) {
        mutProb_.registerCacheableCheckpointTarget(this);
        mutProb_.add_post_change_listener([=, this]() { this->setDirty(); });

        lossProb_.registerCacheableCheckpointTarget(this);
        lossProb_.add_post_change_listener([=, this]() { this->setDirty(); });

        intp_.registerCacheableCheckpointTarget(this);
        intp_.add_set_dirty_listener([=, this]() { this->setDirty(); });

        value_.setZero();
        this->setDirty();
        this->value();
    }



    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    core::datatypes::LogProbabilityTransitionMatrix<2> NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            core::datatypes::SquareMatrix<double, 2> t_mat;
            t_mat <<    1 - mutProb_.value(),   mutProb_.value(),
                    lossProb_.value(),      1 - lossProb_.value();

            auto tmp = t_mat;
            value_ = tmp * intp_.value()(1);

            for (int i = 2; i <= MaxTransmissions; ++i) {
                tmp = tmp * t_mat;
                value_ += tmp * intp_.value()(i);
            }

            value_ = this->value_.array().log();
            this->setClean();
        }

        return value_;

    }


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    double NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLogLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        if(ps.size() > 1) {
            return -std::numeric_limits<double>::infinity();
        }
        double llik = 0.0;
        auto const &childGenotypes = child.latentGenotype();

        auto const childGenotypesIter = childGenotypes.begin();
        for (auto const &parent : ps) {
            auto const &parentGenotypes = parent->latentGenotype();
            auto const parentGenotypesIter = parentGenotypes.begin();
            for (size_t i = 0; i < parentGenotypes.size(); ++i) {
                // assume loci are ordered the same
                auto const &parentGenotypeAtLocus = (*(parentGenotypesIter + i)).second;
                auto const &childGenotypeAtLocus = (*(childGenotypesIter + i)).second;
                const unsigned int t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                const unsigned int t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                const unsigned int t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                const unsigned int t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());

                // no mutation
                llik += t00 * value()(0,0);

                // mutation
                llik += t01 * value()(0,1);

                // loss
                llik += t10 * value()(1,0);

                // no loss
                llik += t11 * value()(1,1);
            }

//        for (auto const& [locus, parentGenotypeAtLocus] : parentGenotypes) {
//            auto const &childGenotypeAtLocus = childGenotypes.at(locus); // this is bad, slow lookup in tight loop -> could represent locus genetic data as one big vector instead?
//            const unsigned int t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
//            const unsigned int t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
//            const unsigned int t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
//            const unsigned int t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
//
//            // no mutation
//            llik += t00 * value()(0,0);
//
//            // mutation
//            llik += t01 * value()(0,1);
//
//            // loss
//            llik += t10 * value()(1,0);
//
//            // no loss
//            llik += t11 * value()(1,1);
//        }
        }
        return llik;
    }

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    double NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLogLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
    //    assert(ps.size() <= 1);
        if(ps.size() > 1) {
            return -std::numeric_limits<double>::infinity();
        }
        double llik = 0.0;
        auto const &childGenotype = child.latentGenotype();
        for (auto const &parent : ps) {
            auto const &parentGenotypes = parent->latentGenotype();
            for (auto const& [locus, parentGenotypeAtLocus] : parentGenotypes) {
                if (childGenotype.contains(locus)) {
                    auto const &childGenotypeAtLocus = childGenotype.at(locus);
                    const unsigned int t00 = GeneticsImpl::trueNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                    const unsigned int t01 = GeneticsImpl::falsePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                    const unsigned int t10 = GeneticsImpl::falseNegativeCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                    const unsigned int t11 = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());

                    // no mutation
                    llik += t00 * peek()(0,0);

                    // mutation
                    llik += t01 * peek()(0,1);

                    // loss
                    llik += t10 * peek()(1,0);

                    // no loss
                    llik += t11 * peek()(1,1);

                }
            }
        }
        return llik;
    }


    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    double NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::calculateLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        double llik = calculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<double>::infinity() ? exp(llik) : 0;
    }

    template<int MaxTransmissions, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    double NoSuperInfectionMutation<MaxTransmissions, InterTransmissionProbImpl>::peekCalculateLikelihood(const core::containers::Infection<GeneticsImpl> &child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) {
        double llik = peekCalculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<double>::infinity() ? exp(llik) : 0;
    }

}


#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONMUTATION_H
