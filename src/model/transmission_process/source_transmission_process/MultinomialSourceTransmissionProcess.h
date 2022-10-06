//
// Created by Maxwell Murphy on 7/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H


#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Locus.h"
#include "core/io/serialize.h"
#include "core/utils/ProbAnyMissing.h"
#include "core/utils/numerics.h"

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include <cmath>
#include <functional>
#include <utility>


namespace transmission_nets::model::transmission_process {
    using Likelihood = core::computation::Likelihood;
    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    class MultinomialSourceTransmissionProcess : public core::computation::PartialLikelihood {

    public:
        MultinomialSourceTransmissionProcess(std::shared_ptr<COIProbabilityImpl> coiProb,
                                             std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer,
                                             std::vector<std::shared_ptr<core::containers::Locus>>  loci,
                                             const GenotypeParameterMap& genetics);

        Likelihood value() override;
        Likelihood validate();
        std::string identifier() override;

    private:
        friend class core::abstract::Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>>;
        friend class core::abstract::Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>, double>;

        void calculateLocusLogLikelihood(std::shared_ptr<core::containers::Locus> locus);
        void postSaveState();
        void postAcceptState();
        void postRestoreState();

        std::shared_ptr<COIProbabilityImpl> coiProb_;
        std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer_;
        std::vector<std::shared_ptr<core::containers::Locus>> loci_{};
        GenotypeParameterMap genetics_;

        // update alleleFrequencies -> update estimates at locus
        // update founder -> update estimates at locus
        // update COI -> update estimates across loci

        boost::container::flat_map<std::shared_ptr<core::containers::Locus>, int> locusIdxMap_{};
        boost::container::flat_set<std::shared_ptr<core::containers::Locus>> dirtyLoci_{};


        // loci independent conditional on COI
        int totalLoci_;

        // buffers for calculations
        std::vector<Likelihood> coiPartialLlik_{};
        std::vector<Likelihood> prVec_{};

        // stateful -- must be cached
        std::vector<Likelihood> locusLlikBuffer_{};
        std::vector<Likelihood> llikMatrix_{};
        std::vector<std::vector<Likelihood>> locusLlikBufferCache_{};
        std::vector<std::vector<Likelihood>> llikMatrixCache_{};

        std::vector<Likelihood> tmpCalculationVec_{};

        core::utils::probAnyMissingFunctor probAnyMissing_;
    };


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>::MultinomialSourceTransmissionProcess(std::shared_ptr<COIProbabilityImpl> coiProb, std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer, std::vector<std::shared_ptr<core::containers::Locus>>  loci, const GenotypeParameterMap& genetics)
        : coiProb_(std::move(coiProb)), alleleFrequenciesContainer_(std::move(alleleFrequenciesContainer)), loci_(std::move(loci)), genetics_(genetics) {
        value_     = 0;
        totalLoci_ = alleleFrequenciesContainer_->totalLoci();

        llikMatrix_.resize((MAX_COI + 1) * totalLoci_);
        coiPartialLlik_.resize(MAX_COI + 1);
        locusLlikBuffer_.resize(MAX_COI + 1);

        coiProb_->registerCacheableCheckpointTarget(this);
        coiProb_->add_set_dirty_listener([=, this]() {
            this->setDirty();
            this->dirtyLoci_.insert(loci_.begin(), loci_.end());
        });

        int idx = 0;
        for (const auto& locus : loci_) {
            this->dirtyLoci_.insert(locus);
            locusIdxMap_[locus] = idx;
            idx++;

            alleleFrequenciesContainer_->alleleFrequencies(locus)->registerCacheableCheckpointTarget(this);
            alleleFrequenciesContainer_->alleleFrequencies(locus)->add_post_change_listener([=, this]() {
                this->setDirty();
                this->dirtyLoci_.insert(locus);
            });

            genetics_.at(locus)->registerCacheableCheckpointTarget(this);
            genetics_.at(locus)->add_post_change_listener([=, this]() {
                this->setDirty();
                this->dirtyLoci_.insert(locus);
            });
        }

        this->addPostSaveHook([=, this]() { this->postSaveState(); });
        this->addPostRestoreHook([=, this]() { this->postRestoreState(); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        this->setDirty();
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    Likelihood MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>::value() {
        if (this->isDirty()) {
            for (const auto& locus : dirtyLoci_) {
                this->calculateLocusLogLikelihood(locus);
                std::copy(locusLlikBuffer_.begin(), locusLlikBuffer_.end(), llikMatrix_.begin() + locusIdxMap_[locus] * (MAX_COI + 1));
            }

            dirtyLoci_.clear();
            tmpCalculationVec_.clear();

            std::fill(coiPartialLlik_.begin(), coiPartialLlik_.end(), 0.0);
            for (int k = 0; k < MAX_COI + 1; ++k) {
                coiPartialLlik_.at(k) += coiProb_->value()[k];
            }

            for (int j = 0; j < totalLoci_; j++) {
                for (int i = 0; i < MAX_COI + 1; i++) {
                    coiPartialLlik_.at(i) += llikMatrix_.at(j * (MAX_COI + 1) + i);
                }
            }

            for (int l = 0; l < MAX_COI + 1; ++l) {
                tmpCalculationVec_.push_back(coiPartialLlik_.at(l));
            }

//            this->value_ = validate();
            this->value_ = core::utils::logSumExp(tmpCalculationVec_);

            if (std::isnan(this->value_) or this->value_ <= -std::numeric_limits<double>::infinity()) {
                fmt::print("NAN in MultinomialSourceTransmissionProcess::value()\n");
                fmt::print("\ttmpCalculationVec_ = {}\n", core::io::serialize(tmpCalculationVec_));
                fmt::print("Allele frequencies:\n");
                for (const auto& locus : loci_) {
                    fmt::print("\tlocus {} = {}\n", locus->label, core::io::serialize(alleleFrequenciesContainer_->alleleFrequencies(locus)->value()));
                    fmt::print("\tgenotype = {}\n", core::io::serialize(genetics_.at(locus)->value()));
                }
                this->value_ = -std::numeric_limits<double>::infinity();
            }

            #ifdef DEBUG_LIKELIHOOD
                        auto tmp = validate();
                        if (std::abs(tmp - this->value_) > 1) {
                            fmt::print(stderr, "Likelihood mismatch MSTP: {}, {}\n", tmp, this->value_);
                        }
            #endif

            this->setClean();
        }


        return this->value_;
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    Likelihood MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::validate() {
        for (const auto& [locus, idx] : locusIdxMap_) {
            this->calculateLocusLogLikelihood(locus);
            std::copy(locusLlikBuffer_.begin(), locusLlikBuffer_.end(), llikMatrix_.begin() + (idx * (MAX_COI + 1)));
        }

        //        dirtyLoci_.clear();
        tmpCalculationVec_.clear();

        std::fill(coiPartialLlik_.begin(), coiPartialLlik_.end(), 0.0);
        for (int k = 0; k < MAX_COI + 1; ++k) {
            coiPartialLlik_.at(k) += coiProb_->value()[k];
        }

        for (int j = 0; j < totalLoci_; j++) {
            for (int i = 0; i < MAX_COI + 1; i++) {
                coiPartialLlik_.at(i) += llikMatrix_.at(j * (MAX_COI + 1) + i);
            }
        }

        for (int l = 0; l < MAX_COI + 1; ++l) {
            tmpCalculationVec_.push_back(coiPartialLlik_.at(l));
        }

        auto value = core::utils::logSumExp(tmpCalculationVec_);
        return value;
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    std::string MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::identifier() {
        return "MultinomialSourceTransmissionProcess";
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    __attribute__((flatten)) void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLocusLogLikelihood(const std::shared_ptr<core::containers::Locus> locus) {
        const auto& alleleFreqs   = alleleFrequenciesContainer_->alleleFrequencies(locus)->value();
        const auto& genotype      = genetics_.at(locus)->value();
        double constrainedSetProb = 0.0;
        bool zeroProbEvent        = false;

        prVec_.clear();
        prVec_.reserve(alleleFreqs.totalElements());
        for (size_t j = 0; j < alleleFreqs.totalElements(); ++j) {
            if (genotype.allele(j)) {
                prVec_.push_back(alleleFreqs.frequencies(j));
                constrainedSetProb += alleleFreqs.frequencies(j);
                zeroProbEvent = std::abs(alleleFreqs.frequencies(j)) < 1e-12;
            }
        }

        if (constrainedSetProb > 0 and !zeroProbEvent) {
            // Normalize the probability density
            for (auto& k : prVec_) {
                k = k / constrainedSetProb;
            }

            for (unsigned int coi = 0; coi < MAX_COI + 1; coi++) {
                // Prob that after `coi` draws 1 or more alleles are not drawn
                Likelihood pam = probAnyMissing_(prVec_, coi);


                // prob that after `coi` draws all alleles are drawn at least once conditional on all draws come from the constrained set.
                if (pam >= 1) {
                    locusLlikBuffer_[coi] = -std::numeric_limits<Likelihood>::infinity();
                } else {
                    locusLlikBuffer_[coi] = std::log(1 - pam) + std::log(constrainedSetProb) * coi;
                }

            }

        } else {
            std::fill(locusLlikBuffer_.begin(), locusLlikBuffer_.end(), -std::numeric_limits<Likelihood>::infinity());
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
        locusLlikBuffer_ = std::move(locusLlikBufferCache_.back());
        locusLlikBufferCache_.pop_back();
        llikMatrix_ = std::move(llikMatrixCache_.back());
        llikMatrixCache_.pop_back();
    }
}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
