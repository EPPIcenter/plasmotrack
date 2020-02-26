//
// Created by Maxwell Murphy on 2/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONTRANSMISSIONPROCESS_H

#include <core/containers/ParentSet.h>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/Matrix.h"
#include "core/datatypes/Alleles.h"

#include "core/containers/ParentSet.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/computation/transmission_process/ZTMultiplicativeBinomial.h"

template<int MAX_TRANSMISSIONS, int MAX_COI>
class NoSuperInfectionTransmissionProcess : public Computation<LogProbabilityMatrix<MAX_COI>>,
                                            public Observable<NoSuperInfectionTransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>,
                                            public Cacheable<NoSuperInfectionTransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>,
                                            public Checkpointable<NoSuperInfectionTransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>, LogProbabilityMatrix<MAX_COI>> {

//  TODO: Generation Probability -- Currently assumes flat distribution over generation times,
//      which is bad when considering potentially many generations between nodes. Maybe implement
//      truncated distributions over generations
public:
    explicit NoSuperInfectionTransmissionProcess(ZTMultiplicativeBinomial<MAX_COI> &ztmb) : ztmb_(ztmb) {
        ztmb_.add_set_dirty_listener([&]() {
            this->setDirty();
        });
        ztmb_.registerCheckpointTarget(*this);
    };

    LogProbabilityMatrix<MAX_COI> value() noexcept override {
        if (this->isDirty()) {
            this->value_ = ztmb_.value();
            auto tmp = ztmb_.value();

            for (int i = 1; i < MAX_TRANSMISSIONS; ++i) {
                tmp = tmp * ztmb_.value();
                this->value_ += tmp;
            }

            this->value_ = (this->value_ / MAX_TRANSMISSIONS).array().log();
            this->setClean();
        }
        return this->value_;
    };


    template <typename GeneticsImpl>
    double calculateLogLikelihood(Infection<GeneticsImpl>& child, ParentSet<Infection<GeneticsImpl>>& ps) {
        assert(ps.size() == 1);
        double llik = 0.0;
        auto const& childGenotype = child.latentGenotype();
        for (auto const& parent : ps) {
            auto const& parentGenotypes = parent->latentGenotype();
            for (auto const& [locus, parentGenotypeAtLocus] : parentGenotypes) {
                if(childGenotype.contains(locus)) {
                    auto const& childGenotypeAtLocus = childGenotype.at(locus);
                    unsigned int parentAlleleCount = parentGenotypeAtLocus.value().totalPositiveCount();
                    unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                    llik += value()(parentAlleleCount, retainedAlleleCount);
                }
            }
        }

        return llik;
    };

private:
    friend class Checkpointable<NoSuperInfectionTransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>, ProbabilityMatrix<MAX_COI>>;
    friend class Cacheable<NoSuperInfectionTransmissionProcess<MAX_TRANSMISSIONS, MAX_COI>>;

    ZTMultiplicativeBinomial<MAX_COI> &ztmb_;
};


#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONTRANSMISSIONPROCESS_H
