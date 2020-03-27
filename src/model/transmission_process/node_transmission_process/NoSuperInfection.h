//
// Created by Maxwell Murphy on 2/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/Matrix.h"
#include "core/datatypes/Alleles.h"

#include "core/containers/ParentSet.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"


template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
class NoSuperInfection : public Computation<LogProbabilityMatrix<MAX_COI + 1>>,
                         public Observable<NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                         public Cacheable<NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                         public Checkpointable<NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>, LogProbabilityMatrix<
                                 MAX_COI + 1>> {

public:
    explicit NoSuperInfection(COITransitionProbImpl &coitp, InterTransmissionProbImpl &intp);

    LogProbabilityMatrix<MAX_COI + 1> value() noexcept override;

    template<typename GeneticsImpl>
    double calculateLogLikelihood(Infection<GeneticsImpl> &child, ParentSet<Infection<GeneticsImpl>> &ps);

    template<typename GeneticsImpl>
    double peekCalculateLogLikelihood(Infection<GeneticsImpl> &child, ParentSet<Infection<GeneticsImpl>> &ps);

    template<typename GeneticsImpl>
    double calculateLikelihood(Infection<GeneticsImpl> &child, ParentSet<Infection<GeneticsImpl>> &ps) {
        return exp(calculateLogLikelihood(child, ps));
    }

    template<typename GeneticsImpl>
    double peekCalculateLikelihood(Infection<GeneticsImpl> &child, ParentSet<Infection<GeneticsImpl>> &ps) {
        return exp(peekCalculateLogLikelihood(child, ps));
    }

private:
    friend class Checkpointable<NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>, LogProbabilityMatrix<
            MAX_COI + 1>>;

    friend class Cacheable<NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>;

    COITransitionProbImpl &coitp_;
    InterTransmissionProbImpl &intp_;
};

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::NoSuperInfection(
        COITransitionProbImpl &coitp, InterTransmissionProbImpl &intp)
        : coitp_(coitp), intp_(intp) {
    coitp_.registerCacheableCheckpointTarget(*this);
    coitp_.add_set_dirty_listener([&]() { this->setDirty(); });

    intp_.registerCacheableCheckpointTarget(*this);
    intp_.add_set_dirty_listener([&]() { this->setDirty(); });
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
LogProbabilityMatrix<MAX_COI + 1>
NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::value() noexcept {
    if (this->isDirty()) {
        this->value_ = coitp_.value() * intp_.value()(1);
        auto tmp = coitp_.value();

        for (int i = 2; i <= MAX_TRANSMISSIONS; ++i) {
            tmp = tmp * coitp_.value() * intp_.value()(i);
            this->value_ += tmp;
        }

        this->value_ = this->value_.array().log();
        this->setClean();
    }
    return this->value_;
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
template<typename GeneticsImpl>
double
NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::calculateLogLikelihood(
        Infection <GeneticsImpl> &child, ParentSet <Infection<GeneticsImpl>> &ps) {
    assert(ps.size() == 1);
    double llik = 0.0;
    auto const &childGenotype = child.latentGenotype();
    for (auto const &parent : ps) {
        auto const &parentGenotypes = parent->latentGenotype();
        for (auto const& [locus, parentGenotypeAtLocus] : parentGenotypes) {
            if (childGenotype.contains(locus)) {
                auto const &childGenotypeAtLocus = childGenotype.at(locus);
                const unsigned int parentAlleleCount = parentGenotypeAtLocus.value().totalPositiveCount();
                const unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(
                        parentGenotypeAtLocus.value(), childGenotypeAtLocus.value());
                llik += value()(parentAlleleCount, retainedAlleleCount);
            }
        }
    }

    return llik;
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
template<typename GeneticsImpl>
double
NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::peekCalculateLogLikelihood(
        Infection <GeneticsImpl> &child, ParentSet <Infection<GeneticsImpl>> &ps) {
    assert(ps.size() == 1);
    double llik = 0.0;
    auto const &childGenotype = child.latentGenotype();
    for (auto const &parent : ps) {
        auto const &parentGenotypes = parent->latentGenotype();
        for (auto const&[locus, parentGenotypeAtLocus] : parentGenotypes) {
            if (childGenotype.contains(locus)) {
                auto const &childGenotypeAtLocus = childGenotype.at(locus);
                unsigned int parentAlleleCount = parentGenotypeAtLocus.value().totalPositiveCount();
                unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(parentGenotypeAtLocus.value(),
                                                                                   childGenotypeAtLocus.value());
                llik += this->peek()(parentAlleleCount, retainedAlleleCount);
            }
        }
    }

    return llik;
}

#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H
