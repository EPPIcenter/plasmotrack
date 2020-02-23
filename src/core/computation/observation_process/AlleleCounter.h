//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELECOUNTER_H
#define TRANSMISSION_NETWORKS_APP_ALLELECOUNTER_H

#include "core/computation/observation_process/AlleleCounts.h"
#include "core/computation/Computation.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"


template<typename GeneticsImpl>
class AlleleCounter : public Computation<AlleleCounts>,
                      public Observable<AlleleCounter<GeneticsImpl>>,
                      public Cacheable<AlleleCounter<GeneticsImpl>>,
                      public Checkpointable<AlleleCounter<GeneticsImpl>, AlleleCounts> {

    static constexpr auto truePositiveCount = &GeneticsImpl::truePositiveCount;
    static constexpr auto falsePositiveCount = &GeneticsImpl::falsePositiveCount;
    static constexpr auto trueNegativeCount = &GeneticsImpl::trueNegativeCount;
    static constexpr auto falseNegativeCount = &GeneticsImpl::falseNegativeCount;

public:
    AlleleCounts value() noexcept override {
        if (this->isDirty()) {
            value_.true_positive_count = truePositiveCount(latent_genetics_.value(), observed_genetics_.value());
            value_.false_positive_count = falsePositiveCount(latent_genetics_.value(), observed_genetics_.value());
            value_.true_negative_count = trueNegativeCount(latent_genetics_.value(), observed_genetics_.value());
            value_.false_negative_count = falseNegativeCount(latent_genetics_.value(), observed_genetics_.value());
            this->setClean();
        }
        return value_;
    }

    AlleleCounter(Parameter<GeneticsImpl> &latentGenetics, Data<GeneticsImpl> &observedGenetics) : latent_genetics_(
            latentGenetics), observed_genetics_(observedGenetics) {
        latent_genetics_.add_post_change_listener([&]() {
            this->setDirty();
        });
        latent_genetics_.registerCheckpointTarget(*this);
    }

private:
    friend class Checkpointable<AlleleCounter<GeneticsImpl>, AlleleCounts>;
    friend class Cacheable<AlleleCounter<GeneticsImpl>>;

    Parameter<GeneticsImpl>& latent_genetics_;
    Data<GeneticsImpl>& observed_genetics_;

};

#endif //TRANSMISSION_NETWORKS_APP_ALLELECOUNTER_H
