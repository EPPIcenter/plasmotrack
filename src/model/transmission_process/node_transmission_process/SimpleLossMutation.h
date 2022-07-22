//
// Created by Maxwell Murphy on 7/11/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H
#define TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"
#include "model/transmission_process/NetworkBasedTransmissionProcess.h"

#include <memory>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    /*
     * Functional Node -- implements calculateLogLikelihood(child, parent_set)
     * Transmission process is a function of the loss and mutation rate per transmission event.
     * Internally calculates a 2x2 transition matrix M representing the probability of a gain or loss of a single allele, which
     * is then integrated over MAX_TRANSMISSIONS to get the probability of an allele being lost or gained
     */
    class SimpleLossMutation : public core::computation::Computation<std::array<long double, MAX_TRANSMISSIONS * 4>>,
                       public core::abstract::Observable<SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>,
                       public core::abstract::Cacheable<SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>,
                       public core::abstract::Checkpointable<SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>, core::computation::Computation<std::array<long double, MAX_TRANSMISSIONS * 4>>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;

    public:
        explicit SimpleLossMutation(p_ParameterDouble loss_prob, p_ParameterDouble mutation_rate, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb);

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) noexcept;


    private:
        friend class core::abstract::Checkpointable<SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>, core::computation::Computation<std::array<long double, MAX_TRANSMISSIONS * 4>>>;
        friend class core::abstract::Cacheable<SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>>;

        std::array<long double, (MAX_TRANSMISSIONS + 1) * 4> value() noexcept override;

        p_ParameterDouble mutProb_;
        p_ParameterDouble lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> interTransProb_;
    };

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::SimpleLossMutation(p_ParameterDouble loss_prob, p_ParameterDouble mutation_rate, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb) :
            mutProb_(std::move(mutation_rate)),
            lossProb_(std::move(loss_prob)),
            interTransProb_(std::move(interTransmissionProb)) {

        mutProb_->template registerCacheableCheckpointTarget(this);
        mutProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        lossProb_->template registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        interTransProb_->template registerCacheableCheckpointTarget(this);
        interTransProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        this->value_.fill(0.0);

        this->setDirty();
        this->value();
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    std::array<long double, (MAX_TRANSMISSIONS + 1) * 4> SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::value() noexcept {
        if (this->isDirty()) {
            this->value_.fill(0.0);

            /*
             * Sketch of transmission matrix and the indices
             *         0   1
             *       |-------
             *     0 | 0   2
             *     1 | 1   3
             *
             */

            // 0 -> 0
            this->value_[0] = 1.0 - this->mutProb_->value();
            // 1 -> 0
            this->value_[1] = this->lossProb_->value();
            // 0 -> 1
            this->value_[2] = this->mutProb_->value();
            // 1 -> 1
            this->value_[3] = 1.0 - this->lossProb_->value();

            for (int i = 1; i <= MAX_TRANSMISSIONS; i++) {
                int j = i - 1;
                // 0 -> 0
                this->value_[i * 4] = this->value_[j * 4] * this->value_[0] + this->value_[j * 4 + 2] * this->value_[1];
                // 1 -> 0
                this->value_[i * 4 + 1] = this->value_[j * 4 + 1] * this->value_[0] + this->value_[j * 4 + 3] * this->value_[1];
                // 0 -> 1
                this->value_[i * 4 + 2] = this->value_[j * 4] * this->value_[2] + this->value_[j * 4 + 2] * this->value_[3];
                // 1 -> 1
                this->value_[i * 4 + 3] = this->value_[j * 4 + 1] * this->value_[2] + this->value_[j * 4 + 3] * this->value_[3];
            }
            this->setClean();
        }
        return this->value_;
    }

    template<int MAX_TRANSMISSIONS, typename InterTransmissionProbImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> child, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>> &ps) noexcept {
        /*
         * P(child = 0 | parents) = P(parent_1 lost allele)*...*P(parent_n lost allele) -- but this doesn't apply when looking at multiple
         * alleles. Need to consider joint genetics of each parent.
         *
         */

    }

}



#endif//TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H
