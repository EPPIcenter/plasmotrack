//
// Created by Maxwell Murphy on 2/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESS_H

#include <core/datatypes/ParentSet.h>
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/Matrix.h"
#include "core/datatypes/ParentSet.h"
#include "core/datatypes/Infection.h"
#include "core/datatypes/Alleles.h"

#include "core/computation/transmission_process/ZTMultiplicativeBinomial.h"

template<int MAX_TRANSMISSIONS, int MAX_COI>
class TransmissionProcess : public Computation<ProbabilityMatrix<MAX_COI>>,
                            public Observable<TransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>,
                            public Cacheable<TransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>,
                            public Checkpointable<TransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>, ProbabilityMatrix<MAX_COI>> {

//  TODO: Generation Probability -- Currently assumes flat distribution over generation times,
//      which is bad when considering potentially many generations between nodes. Maybe implement
//      truncated distributions over generations
public:
    explicit TransmissionProcess(ZTMultiplicativeBinomial<MAX_COI> &ztmb) : ztmb_(ztmb) {
        ztmb_.add_set_dirty_listener([&]() {
            this->setDirty();
        });
        ztmb_.registerCheckpointTarget(*this);
    };

    ProbabilityMatrix<MAX_COI> value() noexcept override {
        if (this->isDirty()) {
            this->value_ = ztmb_.value();
            auto tmp = ztmb_.value();

            for (int i = 1; i < MAX_TRANSMISSIONS; ++i) {
                tmp = tmp * ztmb_.value();
                this->value_ += tmp;
            }
            this->value_ = this->value_ / MAX_TRANSMISSIONS;
            this->setClean();
        }
        return this->value_;
    };

    template <typename GeneticsImpl>
    double calculateLikelihood(Infection<GeneticsImpl> child, ParentSet<Infection<GeneticsImpl>> ps) {
        int psAlleleCount = ps[0];

        return 0;
    };

private:
    friend class Checkpointable<TransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>, ProbabilityMatrix<MAX_COI>>;

    friend class Cacheable<TransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>;

    ZTMultiplicativeBinomial<MAX_COI> &ztmb_;
};


#endif //TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESS_H
