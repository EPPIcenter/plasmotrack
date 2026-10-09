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
#include "core/containers/UidLess.h"
#include "core/io/serialize.h"
#include "core/utils/ProbAnyMissing.h"
#include "core/utils/numerics.h"

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <utility>


namespace transmission_nets::model::transmission_process {
    using Likelihood = core::computation::Likelihood;
    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    class MultinomialSourceTransmissionProcess : public core::computation::PartialLikelihood {

    public:
        MultinomialSourceTransmissionProcess(std::shared_ptr<COIProbabilityImpl> coiProb,
                                             std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer,
                                             std::vector<std::shared_ptr<core::containers::Locus>> loci,
                                             const GenotypeParameterMap& genetics,
                                             bool null_model = false);

        Likelihood value() override;
        Likelihood peek() noexcept override;
        // The log-likelihood computed from scratch, without reading or changing what value() keeps.
        Likelihood validate();
        std::string identifier() override;

    private:
        friend class Cacheable<MultinomialSourceTransmissionProcess>;
        friend class Checkpointable<MultinomialSourceTransmissionProcess, double>;

        // Into locusLlikBuffer_. With `shared`, through the rows shared by the source terms on this thread.
        void calculateLocusLogLikelihood(std::shared_ptr<core::containers::Locus> locus, bool shared = true);
        void postSaveState(int savedStateId);
        void postAcceptState();
        void postRestoreState(int savedStateId);

        std::shared_ptr<COIProbabilityImpl> coiProb_;
        std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer_;
        std::vector<std::shared_ptr<core::containers::Locus>> loci_{};
        GenotypeParameterMap genetics_;

        // update alleleFrequencies -> update estimates at locus
        // update founder -> update estimates at locus
        // update COI -> update estimates across loci

        boost::container::flat_map<std::shared_ptr<core::containers::Locus>, int, core::containers::UidLess> locusIdxMap_{};
        boost::container::flat_set<std::shared_ptr<core::containers::Locus>, core::containers::UidLess> dirtyLoci_{};


        // loci independent conditional on COI
        int totalLoci_;

        // buffers for calculations
        std::vector<Likelihood> coiPartialLlik_{};
        std::vector<Likelihood> prVec_{};

        // Scratch: one locus' log-likelihood for each COI.
        std::vector<Likelihood> locusLlikBuffer_{};
        // Stateful: every locus' row of log-likelihoods by COI.
        std::vector<Likelihood> llikMatrix_{};

        // The rows summed by COI. A changed row is swapped in rather than every row re-summed, so -inf entries are
        // counted rather than added. Swapping accumulates rounding, so the sums are recomputed after MAX_ROW_UPDATES.
        struct RowSums {
            std::array<Likelihood, MAX_COI + 1> finite{};
            std::array<int, MAX_COI + 1> infinite{};
            int updates = 0;
            bool valid  = false;
        };
        static constexpr int MAX_ROW_UPDATES = 1024;
        RowSums rowSums_{};

        // Rows of llikMatrix_ overwritten since each saved state, with their earlier values, and the row sums from
        // before the first of them, so that a restore puts them back without every save copying the whole matrix.
        // Levels are kept and reused across saves.
        struct SavedRow {
            int idx;
            std::array<Likelihood, MAX_COI + 1> values;
        };
        struct SavedLevel {
            std::vector<SavedRow> rows{};
            std::optional<RowSums> sums{};
        };
        std::vector<SavedLevel> savedLevels_{};
        size_t numSavedLevels_ = 0;
        // Copies locusLlikBuffer_ into row idx and swaps it into the row sums, recording the row's earlier values
        // for the current saved state.
        void writeRow(int idx);

        core::utils::probAnyMissingFunctor probAnyMissing_;

        // A row depends only on the locus' allele frequencies and its genotype. After a frequency proposal every
        // infection's source term recomputes that locus, mostly for genotypes other infections share, so rows are
        // kept by genotype for each locus' frequencies as last seen, shared by the source terms on a thread.
        using Genotype = std::decay_t<decltype(std::declval<GenotypeParameterMap&>().at(std::declval<std::shared_ptr<core::containers::Locus>>())->value())>;
        struct GenotypeHash {
            std::size_t operator()(const Genotype& genotype) const noexcept {
                return genotype.hash();
            }
        };
        struct LocusRows {
            std::vector<double> frequencies{};
            std::unordered_map<Genotype, std::array<Likelihood, MAX_COI + 1>, GenotypeHash> rows{};
        };
        static std::unordered_map<const core::containers::Locus*, LocusRows>& sharedRows() {
            thread_local std::unordered_map<const core::containers::Locus*, LocusRows> rows;
            return rows;
        }

        bool null_model_;
    };


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>::MultinomialSourceTransmissionProcess(
            std::shared_ptr<COIProbabilityImpl> coiProb,
            std::shared_ptr<AlleleFrequencyContainer> alleleFrequenciesContainer,
            std::vector<std::shared_ptr<core::containers::Locus>> loci,
            const GenotypeParameterMap& genetics,
            const bool null_model) : coiProb_(std::move(coiProb)), alleleFrequenciesContainer_(std::move(alleleFrequenciesContainer)), loci_(std::move(loci)), genetics_(genetics), null_model_(null_model) {
        value_ = 0;
        totalLoci_ = alleleFrequenciesContainer_->totalLoci();

        llikMatrix_.resize((MAX_COI + 1) * totalLoci_);
        coiPartialLlik_.resize(MAX_COI + 1);
        locusLlikBuffer_.resize(MAX_COI + 1);

        // The rows do not depend on the COI distribution, which enters only when they are combined.
        coiProb_->registerCacheableCheckpointTarget(this);
        coiProb_->add_set_dirty_listener([=, this]() {
            this->setDirty();
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

        this->addPostSaveHook([=, this](const auto savedStateId) { this->postSaveState(savedStateId); });
        this->addPostRestoreHook([=, this](const auto savedStateId) { this->postRestoreState(savedStateId); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        this->setDirty();
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    Likelihood MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>::value() {
        if (null_model_) {
            this->value_ = 0;
            this->setClean();
            return this->value_;
        }
        if (this->isDirty()) {
            for (const auto& locus : dirtyLoci_) {
                this->calculateLocusLogLikelihood(locus);
                writeRow(locusIdxMap_[locus]);
            }
            dirtyLoci_.clear();

            if (not rowSums_.valid or rowSums_.updates >= MAX_ROW_UPDATES) {
                rowSums_ = RowSums{.valid = true};
                for (int j = 0; j < totalLoci_; j++) {
                    for (int i = 0; i < MAX_COI + 1; i++) {
                        const Likelihood entry = llikMatrix_[j * (MAX_COI + 1) + i];
                        if (entry == -std::numeric_limits<Likelihood>::infinity()) {
                            ++rowSums_.infinite[i];
                        } else {
                            rowSums_.finite[i] += entry;
                        }
                    }
                }
            }

            const auto& coiProbs = coiProb_->value();
            for (int k = 0; k < MAX_COI + 1; ++k) {
                coiPartialLlik_[k] = rowSums_.infinite[k] > 0 ? -std::numeric_limits<Likelihood>::infinity() : coiProbs[k] + rowSums_.finite[k];
            }

            this->value_ = core::utils::logSumExp(coiPartialLlik_);

            if (std::isnan(this->value_) or this->value_ <= -std::numeric_limits<double>::infinity()) {
                fmt::print("NAN in MultinomialSourceTransmissionProcess::value()\n");
                fmt::print("\tcoiPartialLlik_ = {}\n", core::io::serialize(coiPartialLlik_));
                fmt::print("Allele frequencies:\n");
                for (const auto& locus : loci_) {
                    fmt::print("\tlocus {} = {}\n", locus->label, core::io::serialize(alleleFrequenciesContainer_->alleleFrequencies(locus)->value()));
                    fmt::print("\tgenotype = {}\n", genetics_.at(locus)->value().allelesStr());
                }
                fmt::print("{}\n", core::io::serialize_matrix(llikMatrix_, totalLoci_, MAX_COI + 1));
                this->value_ = -std::numeric_limits<double>::infinity();
                exit(1);
            }

// #ifdef DEBUG_LIKELIHOOD
//             auto tmp = validate();
//             if (std::abs(tmp - this->value_) > 1) {
//                 fmt::print(stderr, "Likelihood mismatch MSTP: {}, {}\n", tmp, this->value_);
//             }
// #endif

            this->setClean();
        }


        return this->value_;
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GenotypeParameterMap, int MAX_COI>
    Likelihood MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GenotypeParameterMap, MAX_COI>::peek() noexcept {
        if (null_model_) {
            return 0;
        }
        return this->value_;
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    Likelihood MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::validate() {
        std::vector<Likelihood> total(MAX_COI + 1);
        for (int k = 0; k < MAX_COI + 1; ++k) {
            total[k] = coiProb_->value()[k];
        }
        for (const auto& locus : loci_) {
            this->calculateLocusLogLikelihood(locus, false);
            for (int k = 0; k < MAX_COI + 1; ++k) {
                total[k] += locusLlikBuffer_[k];
            }
        }
        return core::utils::logSumExp(total);
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    std::string MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::identifier() {
        return "MultinomialSourceTransmissionProcess";
    }


    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    __attribute__((flatten)) void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::calculateLocusLogLikelihood(const std::shared_ptr<core::containers::Locus> locus, const bool shared) {
        const auto& alleleFreqs = alleleFrequenciesContainer_->alleleFrequencies(locus)->value();
        const auto& genotype = genetics_.at(locus)->value();

        std::array<Likelihood, MAX_COI + 1>* row = nullptr;
        if (shared) {
            auto& locusRows      = sharedRows()[locus.get()];
            bool sameFrequencies = locusRows.frequencies.size() == alleleFreqs.totalElements();
            for (size_t j = 0; sameFrequencies and j < alleleFreqs.totalElements(); ++j) {
                sameFrequencies = locusRows.frequencies[j] == alleleFreqs.frequencies(j);
            }
            if (not sameFrequencies) {
                locusRows.frequencies.resize(alleleFreqs.totalElements());
                for (size_t j = 0; j < alleleFreqs.totalElements(); ++j) {
                    locusRows.frequencies[j] = alleleFreqs.frequencies(j);
                }
                locusRows.rows.clear();
            } else if (const auto found = locusRows.rows.find(genotype); found != locusRows.rows.end()) {
                std::ranges::copy(found->second, locusLlikBuffer_.begin());
                return;
            }
            row = &locusRows.rows[genotype];
        }

        double constrainedSetProb = 0.0;

        prVec_.clear();
        prVec_.reserve(alleleFreqs.totalElements());
        for (size_t j = 0; j < alleleFreqs.totalElements(); ++j) {
            if (genotype.allele(j)) {
                prVec_.push_back(alleleFreqs.frequencies(j));
                constrainedSetProb += alleleFreqs.frequencies(j);
            }
        }


        if (constrainedSetProb > 0) {
            // Normalize the probability density
            for (auto& k : prVec_) {
                k = k / constrainedSetProb;
            }

            const double logConstrainedSetProb = std::log(constrainedSetProb);

            for (unsigned int coi = 0; coi <= MAX_COI; coi++) {
                // Prob that after `coi` draws 1 or more alleles are not drawn
                double pam = probAnyMissing_(prVec_, coi);

                // prob that after `coi` draws all alleles are drawn at least once conditional on all draws come from the constrained set.
                if (pam >= 1 and coi >= prVec_.size()) {
                    locusLlikBuffer_[coi] = -std::numeric_limits<Likelihood>::infinity();
                    // fmt::print("PAM > 1 in MultinomialSourceTransmissionProcess::calculateLocusLogLikelihood()\n");
                    // fmt::print("\tprVec_ = {}\n", core::io::serialize(prVec_));
                    // fmt::print("\tcoi = {}\n", coi);
                    // fmt::print("\talleleFreqs = {}\n", core::io::serialize(alleleFreqs));
                    // fmt::print("\tgenotype = {}\n", core::io::serialize(genotype));
                    // fmt::print("\tconstrainedSetProb = {}\n", constrainedSetProb);
                    // fmt::print("\tpam = {}\n", pam);
                    // fmt::print("\tzeroProbEvent = {}\n", zeroProbEvent);
                } else {
                    locusLlikBuffer_[coi] = std::log(1 - pam) + logConstrainedSetProb * static_cast<double>(coi);
                }
            }

        } else {
            std::ranges::fill(locusLlikBuffer_, -std::numeric_limits<Likelihood>::infinity());
        }
        if (row != nullptr) {
            std::ranges::copy(locusLlikBuffer_, row->begin());
        }
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::writeRow(const int idx) {
        const auto row = llikMatrix_.begin() + idx * (MAX_COI + 1);
        if (numSavedLevels_ > 0) {
            auto& level = savedLevels_[numSavedLevels_ - 1];
            if (not level.sums) {
                level.sums = rowSums_;
            }
            auto& saved = level.rows.emplace_back();
            saved.idx   = idx;
            std::copy(row, row + MAX_COI + 1, saved.values.begin());
        }
        if (rowSums_.valid) {
            constexpr Likelihood negInf = -std::numeric_limits<Likelihood>::infinity();
            for (int i = 0; i < MAX_COI + 1; i++) {
                if (row[i] == negInf) {
                    --rowSums_.infinite[i];
                } else {
                    rowSums_.finite[i] -= row[i];
                }
                if (locusLlikBuffer_[i] == negInf) {
                    ++rowSums_.infinite[i];
                } else {
                    rowSums_.finite[i] += locusLlikBuffer_[i];
                }
            }
            ++rowSums_.updates;
        }
        std::ranges::copy(locusLlikBuffer_, row);
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postSaveState([[maybe_unused]] int savedStateId) {
        if (numSavedLevels_ == savedLevels_.size()) {
            savedLevels_.emplace_back();
        }
        savedLevels_[numSavedLevels_].rows.clear();
        savedLevels_[numSavedLevels_].sums.reset();
        ++numSavedLevels_;
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postAcceptState() {
        numSavedLevels_ = 0;
    }

    template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename InfectionEventImpl, int MAX_COI>
    void MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEventImpl, MAX_COI>::postRestoreState([[maybe_unused]] int savedStateId) {
        // Replayed newest first, so a row written twice gets its value from before the first write.
        const auto& level = savedLevels_[numSavedLevels_ - 1];
        for (auto it = level.rows.rbegin(); it != level.rows.rend(); ++it) {
            std::ranges::copy(it->values, llikMatrix_.begin() + it->idx * (MAX_COI + 1));
        }
        if (level.sums) {
            rowSums_ = *level.sums;
        }
        --numSavedLevels_;
    }
}// namespace transmission_nets::model::transmission_process

#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H