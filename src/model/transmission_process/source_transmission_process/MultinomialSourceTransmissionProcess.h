//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include <boost/math/special_functions/gamma.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/computation/Computation.h"

#include "core/utils/CombinationsWithRepetitionsGenerator.h"
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
                                         InfectionEventImpl &founder)
            : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer), founder_(founder) {

        value_ = 0;
        coiProb_.registerCacheableCheckpointTarget(this);
        coiProb_.registerDirtyTarget(this);



        for(const auto& locus : alleleFrequenciesContainer_.loci) {
            alleleFrequenciesContainer_.alleleFrequencies(locus).add_pre_change_listener([=]() {
                this->setDirty();
                this->value_ -= calculateLocusLogLikelihood(locus);
                this->dirty_loci.push_back(locus);
            });
            alleleFrequenciesContainer_.alleleFrequencies(locus).registerCacheableCheckpointTarget(this);
            this->dirty_loci.push_back(locus);
        }

        founder_.registerCacheableCheckpointTarget(this);
        founder_.add_post_change_listener([=]() { this->setDirty(); });
        this->setDirty();

        for (int j = 1; j <= coiProb_.value().size(); ++j) {
            logFactorials.push_back(lgamma(j));
        }
    }

    double value() override {
        double llik = 0.0;
        for (auto &locus : dirty_loci) {
            llik += this->calculateLocusLogLikelihood(locus);
        }
        dirty_loci.clear();

        this->value_ += llik;
        this->setClean();

        return this->value_;
    };


private:
    friend class Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>;
    friend class Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>, double>;

    COIProbabilityImpl &coiProb_;
    AlleleFrequencyContainer &alleleFrequenciesContainer_;
    InfectionEventImpl &founder_;

    std::vector<Locus *> dirty_loci{};

    std::vector<double> logFactorials{};

    double calculateLocusLogLikelihood(Locus *locus) {
        // Does not check if locus is valid
        const auto &founderGenotypeAtLocus = founder_.latentGenotype(locus);
        const double totalAlleles = locus->totalAlleles();
        const auto &alleleFrequencies = alleleFrequenciesContainer_.alleleFrequencies(locus).value();
        const auto &genotype = founderGenotypeAtLocus.value();
        const unsigned int minCOI = genotype.totalPositiveCount();

        if (minCOI > MAX_COI) {
           return std::numeric_limits<double>::min();
        }

        std::vector<int> alleleVector(totalAlleles, 0);
        std::vector<int> positiveIndices{};

        // Identify locations of positive indices and initialize to 1 and calculate probability of all 1 allele
        for (int k = 0; k < totalAlleles; ++k) {
            alleleVector.at(k) = genotype.allele(k);
            if (alleleVector.at(k) == 1) {
                positiveIndices.push_back(k);
            }
        }

        std::vector<double> logResults{};
        logResults.push_back(calculateMultinomialLogLikelihood(alleleFrequencies, alleleVector) + log(coiProb_.value()(minCOI)));
        double maxResult = logResults.back();

        unsigned int totalObsAlleles = positiveIndices.size();
        std::vector<int> tmpAlleleVector(alleleVector);

        for (unsigned int j = 1; j < MAX_COI - minCOI; ++j) {
            std::vector<double> coiLogResults{};
            double coiMaxResult = std::numeric_limits<double>::lowest();

            CombinationsWithRepetitionsGenerator cs(totalObsAlleles, j);
            while (!cs.completed) {
                memset(tmpAlleleVector.data(), 0, sizeof(int)*tmpAlleleVector.size());
                cs.next();
                for (const auto& idx : cs.curr) {
                    tmpAlleleVector.at(positiveIndices.at(idx))++;
                }
                double result = calculateMultinomialLogLikelihood(alleleFrequencies, tmpAlleleVector);
                coiMaxResult = std::max(result, coiMaxResult);
                coiLogResults.push_back(result);
            }
            auto logCOIProb = log(coiProb_.value()(minCOI + j));
            logResults.push_back(logSumExpKnownMax(coiLogResults.begin(), coiLogResults.end(), coiMaxResult) + logCOIProb);
            maxResult = std::max(maxResult, logResults.back());
        }

        double llik = logSumExpKnownMax(logResults.begin(), logResults.end(), maxResult);
        return llik;
    }

    double calculateMultinomialLogLikelihood(const Simplex& simplex, const std::vector<int>& alleleCounts) {
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

    double calculateLogLikelihood() {
        double llik = 0.0;
        const auto &founderLatentGenotype = founder_.latentGenotype();
        for (auto const &pair : founderLatentGenotype) {
            llik += calculateLocusLogLikelihood(pair.first);
        }

        return std::isnan(llik) ? -std::numeric_limits<double>::infinity() : llik;
    }
};

#endif //TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
