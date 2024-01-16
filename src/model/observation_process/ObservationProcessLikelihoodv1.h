//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV1_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV1_H

#include <cmath>

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

namespace transmission_nets::model::observation_process {

    /*
    * Observation process to incorporate error.
    */

    template<typename AlleleCounter>
    class ObservationProcessLikelihoodv1 : public core::computation::PartialLikelihood {

    public:
        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;
        ObservationProcessLikelihoodv1(std::shared_ptr<AlleleCounter> totalAlleles,
                                       p_Parameterfloat falsePositiveRate,
                                       p_Parameterfloat falseNegativeRate);

        core::computation::Likelihood value() override;
        std::string identifier() override;

    private:
        std::shared_ptr<AlleleCounter> total_alleles_;
        p_Parameterfloat false_positive_rate_;
        p_Parameterfloat false_negative_rate_;
    };

    template<typename AlleleCounter>
    ObservationProcessLikelihoodv1<AlleleCounter>::ObservationProcessLikelihoodv1(std::shared_ptr<AlleleCounter> totalAlleles,
                                                                                  p_Parameterfloat falsePositiveRate,
                                                                                  p_Parameterfloat falseNegativeRate) : total_alleles_(std::move(totalAlleles)),
                                                                                                                         false_positive_rate_(std::move(falsePositiveRate)),
                                                                                                                         false_negative_rate_(std::move(falseNegativeRate)) {

        total_alleles_->add_set_dirty_listener([=, this]() {
            this->setDirty();
        });
        total_alleles_->registerCacheableCheckpointTarget(this);

        false_positive_rate_->add_post_change_listener([=, this]() { this->setDirty(); });
        false_positive_rate_->registerCacheableCheckpointTarget(this);

        false_negative_rate_->add_post_change_listener([=, this]() { this->setDirty(); });
        false_negative_rate_->registerCacheableCheckpointTarget(this);

        this->setDirty();
        this->value();
    }

    template<typename AlleleCounter>
    std::string ObservationProcessLikelihoodv1<AlleleCounter>::identifier() {
        return {"ObservationProcessLikelihoodv1"};
    }

    template<typename AlleleCounter>
    core::computation::Likelihood ObservationProcessLikelihoodv1<AlleleCounter>::value() {
        if (this->isDirty()) {
            value_ = total_alleles_->value().true_positive_count * log(1 - false_positive_rate_->value()) +
                     total_alleles_->value().true_negative_count * log(1 - false_negative_rate_->value()) +
                     total_alleles_->value().false_positive_count * log(false_positive_rate_->value()) +
                     total_alleles_->value().false_negative_count * log(false_negative_rate_->value());
            this->setClean();
        }

        assert(value_ < std::numeric_limits<core::computation::Likelihood>::infinity());

        return value_;
    }


}// namespace transmission_nets::model::observation_process


#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV1_H
