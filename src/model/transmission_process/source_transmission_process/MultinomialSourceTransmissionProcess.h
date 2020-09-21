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
#include "core/io/serialize.h"
#include "core/computation/Computation.h"
#include "core/containers/Locus.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"


namespace transmission_nets::model::transmission_process {
    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    class MultinomialSourceTransmissionProcess : public core::computation::Computation<double>,
                                                 public core::abstract::Observable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>,
                                                 public core::abstract::Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>,
                                                 public core::abstract::Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>, double> {

    public:
        MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb,
                                             AlleleFrequencyContainer &alleleFrequenciesContainer,
                                             InfectionEventImpl &founder);

        double value() override;

    private:
        friend class core::abstract::Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>>;
        friend class core::abstract::Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>, double>;

        void calculateLocusLogLikelihood(core::containers::Locus *locus);
        void postSaveState();
        void postAcceptState();
        void postRestoreState();

        COIProbabilityImpl &coiProb_;
        AlleleFrequencyContainer &alleleFrequenciesContainer_;
        InfectionEventImpl &founder_;

        // update alleleFrequencies -> update estimates at locus
        // update founder -> update estimates at locus
        // update COI -> update estimates across loci

        boost::container::flat_map<core::containers::Locus *, int> locusIdxMap_{};

        boost::container::flat_set<core::containers::Locus *> dirtyLoci_{};

        // loci independent conditional on COI
        int totalLoci_;

        // buffers for calculations
        std::vector<double> coiPartialLlik_{};
        std::vector<double> prVec_{};

        // stateful -- must be cached
        std::vector<double> locusLlikBuffer_{};
        std::vector<double> llikMatrix_{};
        std::vector<std::vector<double>> locusLlikBufferCache_{};
        std::vector<std::vector<double>> llikMatrixCache_{};


        std::vector<double> tmpCalculationVec_{};

        core::utils::probAnyMissingFunctor probAnyMissing_;


    };


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb, AlleleFrequencyContainer &alleleFrequenciesContainer, InfectionEventImpl &founder)
            : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer), founder_(founder) {
        value_ = 0;
        totalLoci_ = alleleFrequenciesContainer_.alleleFrequencies().size();
        llikMatrix_.resize((MAX_COI + 1) * totalLoci_);
        coiPartialLlik_.resize(MAX_COI + 1);
        locusLlikBuffer_.resize(MAX_COI + 1);


        coiProb_.registerCacheableCheckpointTarget(this);
        coiProb_.add_set_dirty_listener([=, this]() {
            this->setDirty();
            this->dirtyLoci_.insert(founder.loci().begin(), founder.loci().end());
        });

        int idx = 0;
        for (const auto &locus : founder.loci()) {
            this->dirtyLoci_.insert(locus);
            locusIdxMap_[locus] = idx;
            idx++;

            alleleFrequenciesContainer_.alleleFrequencies(locus).registerCacheableCheckpointTarget(this);
            alleleFrequenciesContainer_.alleleFrequencies(locus).add_post_change_listener([=, this]() {
                this->setDirty();
                this->dirtyLoci_.insert(locus);
            });

            founder_.latentGenotype(locus).registerCacheableCheckpointTarget(this);
            founder_.latentGenotype(locus).add_post_change_listener([=, this]() {
                this->setDirty();
                this->dirtyLoci_.insert(locus);
            });
        }

        this->addPostSaveHook([=, this]() { this->postSaveState(); });
        this->addPostRestoreHook([=, this]() { this->postRestoreState(); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        this->setDirty();

    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    double MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::value() {
        if (this->isDirty()) {

            for (const auto &locus : dirtyLoci_) {
                calculateLocusLogLikelihood(locus);
                std::copy(locusLlikBuffer_.begin(), locusLlikBuffer_.end(), llikMatrix_.begin() + locusIdxMap_[locus] * (MAX_COI + 1));
            }

            dirtyLoci_.clear();
            tmpCalculationVec_.clear();

            std::fill(coiPartialLlik_.begin(), coiPartialLlik_.end(), 0.0);
            for (int k = 0; k < MAX_COI + 1; ++k) {
                coiPartialLlik_.at(k) += coiProb_.value()[k];
            }

            for(int j = 0; j < totalLoci_; j++) {
                for(int i = 0; i < MAX_COI + 1; i++) {
                    coiPartialLlik_.at(i) += llikMatrix_.at(j * (MAX_COI + 1) + i);
                }
            }

            for (int l = 0; l < MAX_COI + 1; ++l) {
                tmpCalculationVec_.push_back(coiPartialLlik_.at(l));
            }

            this->value_ = core::utils::logSumExp(tmpCalculationVec_);

            assert(!std::isnan(this->value_));

            this->setClean();
        }
//        std::cout << "STP: " << value_ << std::endl;
        return this->value_;
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLocusLogLikelihood(core::containers::Locus *locus) {
        const auto &alleleFreqs = alleleFrequenciesContainer_.alleleFrequencies(locus).value();
        const auto &genotype = founder_.latentGenotype(locus).value();

        double denominator = 0.0;

        prVec_.clear();
        prVec_.reserve(alleleFreqs.totalElements());
        for (size_t j = 0; j < alleleFreqs.totalElements(); ++j) {
            if (genotype.allele(j)) {
                prVec_.push_back(alleleFreqs.frequencies(j));
                denominator += alleleFreqs.frequencies(j);
            }
        }

        if(denominator > 0) {
            // Normalize the vector
            for (double& k : prVec_) {
                k = k / denominator;
            }

            for (int i = 0; i < MAX_COI + 1; i++) {
                // Prob that after i draws 1 or more alleles are not drawn
                double pam = probAnyMissing_(prVec_, i);

                // prob that after i draws all alleles are drawn at least once conditional on all draws come from the constrained set.
                locusLlikBuffer_.at(i) = log(1 - pam) + log(denominator) * i;
            }

        } else {
            std::fill(locusLlikBuffer_.begin(), locusLlikBuffer_.end(), -std::numeric_limits<double>::infinity());
        }

    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postSaveState() {
        locusLlikBufferCache_.emplace_back(locusLlikBuffer_);
        llikMatrixCache_.emplace_back(llikMatrix_);
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postAcceptState() {
        locusLlikBufferCache_.clear();
        llikMatrixCache_.clear();
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postRestoreState() {
        locusLlikBuffer_ = locusLlikBufferCache_.back();
        llikMatrix_ = llikMatrixCache_.back();
        locusLlikBufferCache_.pop_back();
        llikMatrixCache_.pop_back();
    }
}




#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
