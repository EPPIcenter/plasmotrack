//
// Created by Maxwell Murphy on 2/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/Matrix.h"

#include "core/containers/ParentSet.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"


template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
class NoSuperInfectionNoMutation : public Computation<LogProbabilityTransitionMatrix<MAX_COI + 1>>,
                                   public Observable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                                   public Cacheable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>,
                                   public Checkpointable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>, LogProbabilityTransitionMatrix<
                                 MAX_COI + 1>> {

public:
    explicit NoSuperInfectionNoMutation(COITransitionProbImpl &coitp, InterTransmissionProbImpl &intp);

    LogProbabilityTransitionMatrix<MAX_COI + 1> value() noexcept override;

    template<typename GeneticsImpl>
    double calculateLogLikelihood(const Infection<GeneticsImpl> &child, const ParentSet<Infection<GeneticsImpl>> &ps);

    template<typename GeneticsImpl>
    double peekCalculateLogLikelihood(const Infection<GeneticsImpl> &child, const ParentSet<Infection<GeneticsImpl>> &ps);

    template<typename GeneticsImpl>
    double calculateLikelihood(const Infection<GeneticsImpl> &child, const ParentSet<Infection<GeneticsImpl>> &ps) {
        double llik = calculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<double>::infinity() ? exp(llik) : 0;
    }

    template<typename GeneticsImpl>
    double peekCalculateLikelihood(const Infection<GeneticsImpl> &child, const ParentSet<Infection<GeneticsImpl>> &ps) {
        double llik = peekCalculateLogLikelihood(child, ps);
        return llik > -std::numeric_limits<double>::infinity() ? exp(llik) : 0;
    }

private:
    friend class Checkpointable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>, LogProbabilityTransitionMatrix<
            MAX_COI + 1>>;

    friend class Cacheable<NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>>;

    COITransitionProbImpl &coitp_;
    InterTransmissionProbImpl &intp_;
};

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::NoSuperInfectionNoMutation(
        COITransitionProbImpl &coitp, InterTransmissionProbImpl &intp)
        : coitp_(coitp), intp_(intp) {
    coitp_.registerCacheableCheckpointTarget(this);
    coitp_.add_set_dirty_listener([=, this]() { this->setDirty(); });

    intp_.registerCacheableCheckpointTarget(this);
    intp_.add_set_dirty_listener([=, this]() { this->setDirty(); });
    this->setDirty();
    this->value();
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
LogProbabilityTransitionMatrix<MAX_COI + 1>
NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::value() noexcept {
    if (this->isDirty()) {
        auto tmp = coitp_.value();
        this->value_ = tmp * intp_.value()(1);

        for (int i = 2; i <= MAX_TRANSMISSIONS; ++i) {
            tmp = tmp * coitp_.value();
            this->value_ += tmp * intp_.value()(i);
        }

        this->value_ = this->value_.array().log();
        this->setClean();
    }
    return this->value_;
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
template<typename GeneticsImpl>
double
NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::calculateLogLikelihood(
        const Infection <GeneticsImpl> &child, const ParentSet <Infection<GeneticsImpl>> &ps) {
    assert(ps.size() <= 1);
    double llik = 0.0;
    auto const &childGenotype = child.latentGenotype();
    for (auto const &parent : ps) {
        auto const &parentGenotypes = parent->latentGenotype();
        for (auto const& [locus, parentGenotypeAtLocus] : parentGenotypes) {
            if (childGenotype.contains(locus)) {
                auto const &childGenotypeAtLocus = childGenotype.at(locus);
                const unsigned int parentAlleleCount = parentGenotypeAtLocus.value().totalPositiveCount();
                const unsigned int retainedAlleleCount = GeneticsImpl::truePositiveCount(
                        parentGenotypeAtLocus.value(), childGenotypeAtLocus.value()
                        );
                llik += value()(parentAlleleCount, retainedAlleleCount);
            }
        }
    }
    return llik;
}

template<int MAX_COI, int MAX_TRANSMISSIONS, typename COITransitionProbImpl, typename InterTransmissionProbImpl>
template<typename GeneticsImpl>
double
NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>::peekCalculateLogLikelihood(
        const Infection <GeneticsImpl> &child, const ParentSet <Infection<GeneticsImpl>> &ps) {
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

#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTIONNOMUTATION_H
