//
// Created by Maxwell Murphy on 7/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include <cmath>
#include <functional>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include "core/utils/ProbAnyMissing.h"
#include "core/utils/numerics.h"
#include "core/computation/Computation.h"
#include "core/containers/Locus.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"


template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
class MultinomialSourceTransmissionProcess : public Computation<double>,
                                             public Observable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>,
                                             public Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>,
                                             public Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>, double> {
    using CallbackType = std::function<void()>;
    CREATE_EVENT(save_state, CallbackType)
    CREATE_EVENT(accept_state, CallbackType)
    CREATE_EVENT(restore_state, CallbackType)

public:
    MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb,
                                         AlleleFrequencyContainer &alleleFrequenciesContainer,
                                         InfectionEventImpl &founder);

    double value() override;

private:
    friend class Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>;
    friend class Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>, double>;

    COIProbabilityImpl &coiProb_;
    AlleleFrequencyContainer &alleleFrequenciesContainer_;
    InfectionEventImpl &founder_;

    // update alleleFrequencies -> update estimates at locus
    // update founder -> update estimates at locus
    // update COI -> update estimates across loci


    boost::container::flat_set<Locus *> dirtyLoci{};

    // loci independent conditional on COI
    std::vector<double> llikMatrix{};
    boost::container::flat_map<Locus *, int> locusIdxMap{};
    int totalLoci;

    std::vector<double> coiPartialLlik{};
    std::vector<double> tmpLocusLlik{};
    std::vector<double> prVec{};

    std::vector<double> tmpCalculationVec{};

    probAnyMissingFunctor probAnyMissing;

    void calculateLocusLogLikelihood(Locus *locus);

};


template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb, AlleleFrequencyContainer &alleleFrequenciesContainer, InfectionEventImpl &founder)
        : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer), founder_(founder) {
    value_ = 0;
    totalLoci = alleleFrequenciesContainer_.alleleFrequencies().size();
    llikMatrix.resize((MAX_COI + 1) * totalLoci);
    coiPartialLlik.resize(MAX_COI + 1);
    tmpLocusLlik.resize(MAX_COI + 1);


    coiProb_.registerCacheableCheckpointTarget(this);
    coiProb_.add_set_dirty_listener([=, this]() {
      this->setDirty();
      this->dirtyLoci.insert(founder.loci().begin(), founder.loci().end());
    });

    int idx = 0;
    for (const auto &locus : founder.loci()) {
        this->dirtyLoci.insert(locus);
        locusIdxMap[locus] = idx;
        idx++;

        alleleFrequenciesContainer_.alleleFrequencies(locus).registerCacheableCheckpointTarget(this);
        alleleFrequenciesContainer_.alleleFrequencies(locus).add_post_change_listener([=, this]() {
          this->setDirty();
          this->dirtyLoci.insert(locus);
        });

        founder_.latentGenotype(locus).registerCacheableCheckpointTarget(this);
        founder_.latentGenotype(locus).add_post_change_listener([=, this]() {
          this->setDirty();
          this->dirtyLoci.insert(locus);
        });
    }

    this->setDirty();

}


template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::value() {
    if (this->isDirty()) {

        for (const auto &locus : dirtyLoci) {
            calculateLocusLogLikelihood(locus);
//            const auto res = calculateLocusLogLikelihood(locus);
            std::copy(tmpLocusLlik.begin(), tmpLocusLlik.end(), llikMatrix.begin() + locusIdxMap[locus] * (MAX_COI + 1));
        }

        dirtyLoci.clear();
        tmpCalculationVec.clear();

        std::fill(coiPartialLlik.begin(), coiPartialLlik.end(), 0.0);
        for (int k = 0; k < MAX_COI + 1; ++k) {
            coiPartialLlik.at(k) += coiProb_.value()(k);
        }

        for(int j = 0; j < totalLoci; j++) {
            for(int i = 0; i < MAX_COI + 1; i++) {
                coiPartialLlik.at(i) += llikMatrix.at(j * (MAX_COI + 1) + i);
            }
        }

        for (int l = 0; l < MAX_COI + 1; ++l) {
            tmpCalculationVec.push_back(coiPartialLlik.at(l));
        }

        this->value_ = logSumExp(tmpCalculationVec);

        this->setClean();
    }

    return this->value_;
}


template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLocusLogLikelihood(Locus *locus) {
    const auto &alleleFreqs = alleleFrequenciesContainer_.alleleFrequencies(locus).value();
    const auto &genotype = founder_.latentGenotype(locus).value();

    double denominator = 0.0;

    prVec.clear();
    prVec.reserve(alleleFreqs.totalElements());
    for (size_t j = 0; j < alleleFreqs.totalElements(); ++j) {
        if (genotype.allele(j)) {
            prVec.push_back(alleleFreqs.frequencies(j));
            denominator += alleleFreqs.frequencies(j);
        }
    }

    for (double & k : prVec) {
        k = k / denominator;
    }

//    std::fill(tmpLocusLlik.begin(), tmpLocusLlik.end(), 0.0);
    if(denominator > 0) {
        for (int i = 0; i < MAX_COI + 1; i++) {
            double pam = probAnyMissing(prVec, i);
            tmpLocusLlik.at(i) = log(1 - pam) + log(denominator) * i;
        }
    } else {
        std::fill(tmpLocusLlik.begin(), tmpLocusLlik.end(), -std::numeric_limits<double>::infinity());
    }

//    return tmpLocusLlik;
}


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
