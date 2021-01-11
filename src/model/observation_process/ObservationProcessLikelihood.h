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

        core::computation::Likelihood value() override;
        std::string identifier() override;

    private:
        AlleleCounter &total_alleles_;
        core::parameters::Parameter<double> &false_positive_rate_;
        core::parameters::Parameter<double> &false_negative_rate_;
    };

    template<typename AlleleCounter>
    ObservationProcessLikelihood<AlleleCounter>::ObservationProcessLikelihood(AlleleCounter &totalAlleles,
                                                                              core::parameters::Parameter<double> &falsePositiveRate,
                                                                              core::parameters::Parameter<double> &falseNegativeRate) : total_alleles_(totalAlleles),
            false_positive_rate_(falsePositiveRate),
            false_negative_rate_(falseNegativeRate) {

        total_alleles_.add_set_dirty_listener([=, this]() {
            this->setDirty();
        });
        total_alleles_.registerCacheableCheckpointTarget(this);

        false_positive_rate_.add_post_change_listener([=, this]() { this->setDirty(); });
        false_positive_rate_.registerCacheableCheckpointTarget(this);

        false_negative_rate_.add_post_change_listener([=, this]() { this->setDirty(); });
        false_negative_rate_.registerCacheableCheckpointTarget(this);

        this->setDirty();
        this->value();
    }

    template<typename AlleleCounter>
    std::string ObservationProcessLikelihood<AlleleCounter>::identifier() {
        return std::string("ObservationProcessLikelihood");
    }

    template<typename AlleleCounter>
    core::computation::Likelihood ObservationProcessLikelihood<AlleleCounter>::value() {
        if (this->isDirty()) {
//            std::cout << "OBS: " << value_ << " | ";
            value_ = total_alleles_.value().true_positive_count * log(1 - false_positive_rate_.value()) +
                     total_alleles_.value().true_negative_count * log(1 - false_negative_rate_.value()) +
                     total_alleles_.value().false_positive_count * log(false_positive_rate_.value()) +
                     total_alleles_.value().false_negative_count * log(false_negative_rate_.value());
//            std::cout << value_ << std::endl;
//            std::cout << total_alleles_.value().true_positive_count << ", ";
//            std::cout << total_alleles_.value().true_negative_count << ", ";
//            std::cout << total_alleles_.value().false_positive_count << ", ";
//            std::cout << total_alleles_.value().false_negative_count << ", ";
//            std::cout << total_alleles_.value().true_positive_count +
//                         total_alleles_.value().true_negative_count +
//                         total_alleles_.value().false_positive_count +
//                         total_alleles_.value().false_negative_count << std::endl;

            this->setClean();
        }

        assert(value_ < std::numeric_limits<core::computation::Likelihood>::infinity());

        return value_;
    }


}


#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOOD_H
