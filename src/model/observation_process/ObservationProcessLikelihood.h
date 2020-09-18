//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H

#include <cmath>

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

namespace transmission_nets::model::observation_process {


    template<typename AlleleCounter>
    class ObservationProcessLikelihood : public core::computation::PartialLikelihood {

    public:
        ObservationProcessLikelihood(AlleleCounter &totalAlleles,
                                     core::parameters::Parameter<double> &falsePositiveRate,
                                     core::parameters::Parameter<double> &falseNegativeRate);

        double value() override;

    private:
        AlleleCounter &total_alleles;
        core::parameters::Parameter<double> &false_positive_rate_;
        core::parameters::Parameter<double> &false_negative_rate_;
    };

    template<typename AlleleCounter>
    ObservationProcessLikelihood<AlleleCounter>::ObservationProcessLikelihood(AlleleCounter &totalAlleles,
                                                                              core::parameters::Parameter<double> &falsePositiveRate,
                                                                              core::parameters::Parameter<double> &falseNegativeRate) :
            total_alleles(totalAlleles),
            false_positive_rate_(falsePositiveRate),
            false_negative_rate_(falseNegativeRate) {
        total_alleles.add_set_dirty_listener([=, this]() { this->setDirty(); });
        total_alleles.registerCacheableCheckpointTarget(this);

        false_positive_rate_.add_post_change_listener([=, this]() { this->setDirty(); });
        false_positive_rate_.registerCacheableCheckpointTarget(this);

        false_negative_rate_.add_post_change_listener([=, this]() { this->setDirty(); });
        false_negative_rate_.registerCacheableCheckpointTarget(this);

        this->setDirty();
        this->value();
    }

    template<typename AlleleCounter>
    double ObservationProcessLikelihood<AlleleCounter>::value() {
        if (this->isDirty()) {
            value_ = total_alleles.value().true_positive_count * log(1 - false_positive_rate_.value()) +
                     total_alleles.value().true_negative_count * log(1 - false_negative_rate_.value()) +
                     total_alleles.value().false_positive_count * log(false_positive_rate_.value()) +
                     total_alleles.value().false_negative_count * log(false_negative_rate_.value());
            this->setClean();
        }

        assert(this->value_ < std::numeric_limits<double>::infinity());

        return value_;
    }


}


#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H
