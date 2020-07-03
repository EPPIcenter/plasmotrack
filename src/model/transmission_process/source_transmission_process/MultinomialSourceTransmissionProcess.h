//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include <functional>

#include <boost/container/flat_set.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"

#include "core/computation/Computation.h"

#include "core/utils/CombinationsWithRepetitionsGenerator.h"
#include "core/utils/io/serialize.h"
#include "core/utils/numerics.h"


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

    CombinationsWithRepetitionsGenerator cs_;

    boost::container::flat_set<Locus *> dirty_loci{};
    boost::container::flat_map<Locus *, double> loci_llik{};
    boost::container::flat_map<Locus *, unsigned int> locusCOI{};

    std::vector<double> logFactorials{};
    std::vector<double> logResults{};
    std::vector<double> coiLogResults{};

    std::vector<int> alleleVector{};
    std::vector<int> tmpAlleleVector{};
    std::vector<int> positiveIndices{};

    double calculateLocusLogLikelihood(Locus *locus);

    double calculateMultinomialLogLikelihood(const Simplex &simplex, const std::vector<int> &alleleCounts);

    double calculateLogLikelihood();

    void updateLocusCOI(Locus *locus);
};

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb, AlleleFrequencyContainer &alleleFrequenciesContainer, InfectionEventImpl &founder)
        : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer), founder_(founder) {

    value_ = 0;
    coiProb_.registerCacheableCheckpointTarget(this);
    coiProb_.add_set_dirty_listener([=]() {
      this->setDirty();
      this->dirty_loci.insert(founder.loci().begin(), founder.loci().end());
    });

    for (const auto &locus : founder.loci()) {
        this->dirty_loci.insert(locus);
        this->loci_llik.emplace(locus, 0);
        this->updateLocusCOI(locus);

        alleleFrequenciesContainer_.alleleFrequencies(locus).registerCacheableCheckpointTarget(this);
        alleleFrequenciesContainer_.alleleFrequencies(locus).add_post_change_listener([=]() {
          this->setDirty();
          this->dirty_loci.insert(locus);
        });

        founder_.latentGenotype(locus).registerCacheableCheckpointTarget(this);
        founder_.latentGenotype(locus).add_post_change_listener([=]() {
          this->setDirty();
          this->dirty_loci.insert(locus);
          this->updateLocusCOI(locus);
        });
    }

    logFactorials.reserve(coiProb_.value().size());
    for (int j = 1; j <= coiProb_.value().size(); ++j) {
        logFactorials.push_back(lgamma(j));
    }
    this->value();
}

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::value() {
    if (this->isDirty()) {
        for (const auto &locus : dirty_loci) {
            loci_llik.at(locus) = calculateLocusLogLikelihood(locus);
        }

        dirty_loci.clear();
        this->value_ = 0.0;

        for (const auto &[locus, val] : loci_llik) {
            this->value_ += val;
        }
        this->setClean();
    }

    //        std::cout << "STP: " << this->value_ << std::endl;
    return this->value_;
}

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLocusLogLikelihood(Locus *locus) {
    // Does not check if locus is valid
    const auto &founderGenotypeAtLocus = founder_.latentGenotype(locus);
    const double totalAlleles = locus->totalAlleles();
    const auto &alleleFrequencies = alleleFrequenciesContainer_.alleleFrequencies(locus).value();
    const auto &genotype = founderGenotypeAtLocus.value();
    const unsigned int minCOI = std::min_element(locusCOI.begin(), locusCOI.end(), [](const auto& l, const auto& r) { return l.second < r.second; })->second;

    if (minCOI > MAX_COI or minCOI == 0) {
        return -std::numeric_limits<double>::infinity();
    }

    alleleVector.resize(totalAlleles);
    memset(alleleVector.data(), 0, sizeof(int) * alleleVector.size());

    positiveIndices.clear();

    // Identify locations of positive indices, initialize to 1, calculate probability
    for (int k = 0; k < totalAlleles; ++k) {
        alleleVector.at(k) = genotype.allele(k);
        if (alleleVector.at(k) == 1) {
            positiveIndices.push_back(k);
        }
    }

    logResults.clear();
    logResults.push_back(calculateMultinomialLogLikelihood(alleleFrequencies, alleleVector) + log(coiProb_.value()(minCOI)));
    double maxResult = logResults.back();

    unsigned int totalObsAlleles = positiveIndices.size();
    tmpAlleleVector = alleleVector;

    for (unsigned int j = 1; j < MAX_COI - minCOI; ++j) {
        coiLogResults.clear();
        double coiMaxResult = std::numeric_limits<double>::lowest();

        cs_.reset(totalObsAlleles, j);
        while (!cs_.completed) {
            tmpAlleleVector = alleleVector;
            cs_.next();
            for (const auto &idx : cs_.curr) {
                tmpAlleleVector.at(positiveIndices.at(idx))++;
            }
            double result = calculateMultinomialLogLikelihood(alleleFrequencies, tmpAlleleVector);
            coiMaxResult = std::max(result, coiMaxResult);
            coiLogResults.push_back(result);
        }
        double logCOIProb = log(coiProb_.value()(minCOI + j));
        logResults.push_back(logSumExpKnownMax(coiLogResults.begin(), coiLogResults.end(), coiMaxResult) + logCOIProb);
        maxResult = std::max(maxResult, logResults.back());
    }
    double llik = logSumExpKnownMax(logResults.begin(), logResults.end(), maxResult);
    return llik;
}

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateMultinomialLogLikelihood(const Simplex &simplex, const std::vector<int> &alleleCounts) {
    int n = 0;
    double llik = 0;
    const double totalElements = simplex.totalElements();
    for (int j = 0; j < totalElements; ++j) {
        n += alleleCounts.at(j);
        llik += log(simplex.frequencies(j)) * alleleCounts.at(j);
        llik -= logFactorials.at(alleleCounts.at(j));
    }
    llik += logFactorials.at(n);
    return llik;
}

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLogLikelihood() {
    double llik = 0.0;
    const auto &founderLatentGenotype = founder_.latentGenotype();
    for (auto const &pair : founderLatentGenotype) {
        llik += calculateLocusLogLikelihood(pair.first);
    }
    return std::isnan(llik) ? -std::numeric_limits<double>::infinity() : llik;
}

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::updateLocusCOI(Locus* locus) {
    locusCOI[locus] = founder_.latentGenotype(locus).value().totalPositiveCount();
}


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
