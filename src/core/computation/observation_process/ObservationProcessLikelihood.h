//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H

#include <cmath>

#include "core/computation/PartialLikelihood.h"
#include "core/parameter/Parameter.h"

template<typename AlleleCounter>
class ObservationProcessLikelihood : public PartialLikelihood {

public:
    ObservationProcessLikelihood(const std::string &id,
                                 AlleleCounter &totalAlleles,
                                 Parameter<double> &falsePositiveRate,
                                 Parameter<double> &falseNegativeRate) :
            PartialLikelihood(id),
            total_alleles(totalAlleles),
            false_positive_rate_(falsePositiveRate),
            false_negative_rate_(falseNegativeRate) {
        total_alleles.add_set_dirty_listener([&]() {this->setDirty();});
        total_alleles.registerCheckpointTarget(*this);

        false_positive_rate_.add_post_change_listener([&]() {this->setDirty();});
        false_positive_rate_.registerCheckpointTarget(*this);

        false_negative_rate_.add_post_change_listener([&]() {this->setDirty();});
        false_negative_rate_.registerCheckpointTarget(*this);

    }

    double value() override {
        if (this->isDirty()) {
            value_ = total_alleles.value().true_positive_count * log(1 - false_positive_rate_.value()) +
                     total_alleles.value().true_negative_count * log(1 - false_negative_rate_.value()) +
                     total_alleles.value().false_positive_count * log(false_positive_rate_.value()) +
                     total_alleles.value().false_negative_count * log(false_negative_rate_.value());
            this->setClean();
        }
        return value_;
    }

private:
    AlleleCounter &total_alleles;
    Parameter<double> &false_positive_rate_;
    Parameter<double> &false_negative_rate_;
};


#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H
