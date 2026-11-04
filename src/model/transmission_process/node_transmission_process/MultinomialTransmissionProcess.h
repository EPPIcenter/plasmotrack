//
// Created by Maxwell Murphy on 11/17/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"

#include "core/utils/ProbAnyMissing.h"
#include "core/utils/numerics.h"

#include <boost/math/special_functions/binomial.hpp>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;
    using Probability = core::computation::Probability;

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    class MultinomialTransmissionProcess : public core::computation::Computation<std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>>,
                                           public core::abstract::Observable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>>,
                                           public core::abstract::Cacheable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>>,
                                           public core::abstract::Checkpointable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
        using p_ParameterArray = std::shared_ptr<core::parameters::Parameter<std::array<Probability, MAX_PARENTS + 1>>>;

        using p_Infection = std::shared_ptr<core::containers::Infection<GeneticsImpl>>;
        using ParentSet   = core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>;

        using p_ParentSetSizePrior = std::shared_ptr<ParentSetSizePriorImpl>;

    public:
        explicit MultinomialTransmissionProcess(p_ParameterDouble mean_strains_transmitted);

        std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)> value() override;

        Likelihood probNumStrains(const int num_strains, const int num_parents) {
            // value() returns the whole table by value; read the cached table instead of copying it per lookup.
            if (this->isDirty()) {
                this->value();
            }
            return this->value_[(num_parents - 1) * MAX_STRAINS + (num_strains - 1)];
        }

        // The term for the observed parents alone, with the parent-set-size prior. The Model evaluates parent sets
        // in batches (below) and sums the latent parent out in LatentParentMarginal; this is the plain form the tests
        // check those against.
        Likelihood calculateLogLikelihood(p_Infection infection, const ParentSet& parentSet, p_ParentSetSizePrior psp);

        // Batch evaluation of many parent sets of one child, with genome identities resolved once rather than per
        // parent set. Call beginBatch(), resolve ids with genomeOf() / locusOrderOf(), then call parentSetSums(); the
        // ids stay valid until the next beginBatch() or calculateLogLikelihood() call, and the genotypes involved must
        // not change in between.
        using GenomeId = std::uint32_t;

        struct ResolvedGenomes {
            GenomeId childLocusOrder = 0;
            GenomeId child           = 0;
            GenomeId latentParent    = 0;
            std::array<GenomeId, MAX_PARENTS> parents{};// in the parent set's iteration order
        };

        void beginBatch() {
            reserveCacheCapacity();
        }
        GenomeId genomeOf(const p_Infection& infection) {
            return genomeId(infection);
        }
        GenomeId locusOrderOf(const p_Infection& infection) {
            return locusOrderId(infection);
        }

        // How each locus' parent allele pool is built: from the observed parents, from them plus a latent parent, or
        // from a latent parent alone.
        enum class Variant : unsigned char {
            ObservedParents,
            ObservedAndLatentParents,
            LatentParentOnly
        };

        using StrainLogLikelihoods = std::array<Likelihood, MAX_STRAINS>;

        // Contribution of one locus to each strain count. `impossible` means the whole parent set has zero likelihood.
        struct LocusTerm {
            bool impossible = false;
            StrainLogLikelihoods terms{};
        };

        using Locus = typename core::containers::Infection<GeneticsImpl>::LocusType;

        // Locus terms summed over a child's loci, kept so that one locus' term can be swapped for another: -inf terms
        // and impossible loci are counted rather than added, so removing them is exact.
        struct StrainSums {
            StrainLogLikelihoods finite{};
            std::array<std::uint16_t, MAX_STRAINS> infinite{};
            std::uint16_t impossible = 0;
            // Set when summing stopped at this locus, the first that makes the parent set impossible. The sums then
            // cover only the loci before it, so no term can be swapped in or out; the parent set stays impossible
            // until this locus' term changes.
            const Locus* stoppedAt = nullptr;

            void add(const LocusTerm& term) {
                if (term.impossible) {
                    ++impossible;
                    return;
                }
                for (unsigned int idx = 0; idx < MAX_STRAINS; ++idx) {
                    if (term.terms[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                        ++infinite[idx];
                    } else {
                        finite[idx] += term.terms[idx];
                    }
                }
            }

            void remove(const LocusTerm& term) {
                if (term.impossible) {
                    --impossible;
                    return;
                }
                for (unsigned int idx = 0; idx < MAX_STRAINS; ++idx) {
                    if (term.terms[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                        --infinite[idx];
                    } else {
                        finite[idx] -= term.terms[idx];
                    }
                }
            }
        };

        // A parent set's locus terms summed over the child's loci. parentSet is null for the latent-parent-only
        // variant and latentParent null for the observed-parents variant. The reference stays valid until the next
        // beginBatch() or calculateLogLikelihood() call.
        const StrainSums& parentSetSums(Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, const ResolvedGenomes& resolved);

        // One locus' term, from the genotypes at that locus: the child's, the observed parents' in parent-set order,
        // and the latent parent's (null for the observed-parents variant). Valid until the next beginBatch() or
        // calculateLogLikelihood() call.
        const LocusTerm& locusTerm(Variant variant, const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents, const GeneticsImpl* latentParent);

        // As locusTerm, computed afresh and not cached: for callers that keep their own cache of the terms they need,
        // whose terms would only push reused ones out of this cache.
        LocusTerm uncachedLocusTerm(Variant variant, const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents, const GeneticsImpl* latentParent);

        // The parent set's log-likelihood marginalised over strain counts, before the source and parent-set-size
        // terms. numParents counts the latent parent.
        Likelihood strainSummedLogLikelihood(const StrainSums& sums, size_t numParents);


    private:
        friend class core::abstract::Checkpointable<MultinomialTransmissionProcess, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>>;
        friend class core::abstract::Cacheable<MultinomialTransmissionProcess>;

        // A locus term depends only on these genotypes, so identical genotypes (common in clonal clusters)
        // share one evaluation. genotypes[0] is the child, then the parents in parent-set order, latent parent last.
        struct LocusTermKey {
            Variant variant;
            unsigned int numParents;
            std::array<GeneticsImpl, MAX_PARENTS + 2> genotypes{};

            bool operator==(const LocusTermKey& rhs) const = default;
        };

        struct LocusTermKeyHash {
            std::size_t operator()(const LocusTermKey& key) const noexcept {
                std::size_t h = static_cast<std::size_t>(key.variant) * 31 + key.numParents;
                for (unsigned int i = 0; i <= key.numParents; ++i) {
                    h = (h * 0x100000001b3ULL) ^ key.genotypes[i].hash();
                }
                return h;
            }
        };

        const LocusTerm& locusTerm(const LocusTermKey& key);

        LocusTerm computeLocusTerm(const LocusTermKey& key);

        // Genome-level memo. Re-evaluating a parent set means summing locus terms over every locus, but in a clonal
        // cluster most parent sets are made of the same few genomes, so the per-locus sums are cached by the identity
        // of the whole genomes involved. Ids are interned from genome content; an infection's id is recomputed only
        // when its latentGenomeVersion() changes.
        struct InfectionGenome {
            std::uint64_t version;
            GenomeId id;
        };

        struct GenomeHash {
            std::size_t operator()(const std::vector<GeneticsImpl>& genome) const noexcept {
                std::size_t h = genome.size();
                for (const auto& genotype : genome) {
                    h = (h * 0x100000001b3ULL) ^ genotype.hash();
                }
                return h;
            }
        };

        struct LocusOrderHash {
            std::size_t operator()(const std::vector<const void*>& order) const noexcept {
                std::size_t h = order.size();
                for (const auto* locus : order) {
                    h = (h * 0x100000001b3ULL) ^ std::hash<const void*>{}(locus);
                }
                return h;
            }
        };

        // Summation over loci follows the child's locus order, so that order is part of the key alongside the genomes.
        struct StrainSumKey {
            Variant variant;
            unsigned int numParents;
            GenomeId childLocusOrder;
            std::array<GenomeId, MAX_PARENTS + 2> genomes{};

            bool operator==(const StrainSumKey& rhs) const = default;
        };

        struct StrainSumKeyHash {
            std::size_t operator()(const StrainSumKey& key) const noexcept {
                std::size_t h = (static_cast<std::size_t>(key.variant) * 31 + key.numParents) * 0x9e3779b97f4a7c15ULL ^ key.childLocusOrder;
                for (unsigned int i = 0; i <= key.numParents; ++i) {
                    h = (h * 0x100000001b3ULL) ^ key.genomes[i];
                }
                return h;
            }
        };

        // Bounds on cache memory. Every cache is cleared together, at the start of a calculation, so ids handed out
        // earlier in the same calculation stay valid.
        static constexpr std::size_t MAX_LOCUS_TERM_CACHE_SIZE = 1 << 17;
        static constexpr std::size_t MAX_STRAIN_SUM_CACHE_SIZE = 1 << 19;
        static constexpr std::size_t MAX_INTERNED_GENOMES      = 1 << 14;

        void reserveCacheCapacity();
        GenomeId genomeId(const p_Infection& infection);
        GenomeId locusOrderId(const p_Infection& infection);

        // A parent set's per-locus terms summed over the child's loci, and its strain-summed log-likelihood for the
        // strain-count table version it was computed with.
        struct StrainSumEntry {
            StrainSums sums{};
            std::uint64_t tableVersion = 0;
            Likelihood summed          = 0;
        };

        // `resolved` may be null, in which case genome ids are looked up here.
        StrainSumEntry& strainSums(Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, size_t numParents, const ResolvedGenomes* resolved);

        // Adds the strain-count prior to the strain sums and log-sum-exps over strain counts.
        // Empty if some locus makes the parent set impossible.
        std::optional<Likelihood> strainSummedLogLikelihood(Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, size_t numParents, const ResolvedGenomes* resolved);

        // Body of calculateLogLikelihood, also callable with pre-resolved genome ids.
        Likelihood logLikelihood(Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, const p_ParentSetSizePrior& psp, const ResolvedGenomes* resolved);

        p_ParameterDouble meanStrainsTransmitted_;
        p_ParameterArray probParentSetSize_;

        core::utils::probAnyMissingFunctor probAnyMissing_;

        // Scratch buffers reused across locus evaluations.
        std::vector<size_t> childAlleles_{};
        std::vector<Probability> parentPopFreqs_{};
        std::vector<Likelihood> prVec_{};
        std::vector<GeneticsImpl> genomeScratch_{};

        // Caches are per instance (one per chain); all keys are pure functions of genotype content.
        std::unordered_map<LocusTermKey, LocusTerm, LocusTermKeyHash> locusTerms_{};
        std::unordered_map<StrainSumKey, StrainSumEntry, StrainSumKeyHash> strainSums_{};

        // Changes whenever the strain-count table (value()) is recomputed or restored. 0 is never current.
        std::uint64_t tableVersion_ = 0;
        std::unordered_map<std::vector<GeneticsImpl>, GenomeId, GenomeHash> genomeIds_{};
        std::unordered_map<std::vector<const void*>, GenomeId, LocusOrderHash> locusOrderIds_{};
        std::unordered_map<int, InfectionGenome> infectionGenomes_{};
        std::unordered_map<int, GenomeId> infectionLocusOrders_{};
#ifndef NDEBUG
        std::vector<const std::vector<GeneticsImpl>*> genomesById_{};
#endif
    };

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::MultinomialTransmissionProcess(p_ParameterDouble mean_strains_transmitted) : meanStrainsTransmitted_(std::move(mean_strains_transmitted)) {
        meanStrainsTransmitted_->registerCacheableCheckpointTarget(this);
        meanStrainsTransmitted_->add_post_change_listener([this]() {
            this->setDirty();
        });
        // A restore brings back an earlier table without recomputing it.
        this->addPostRestoreHook([this]([[maybe_unused]] int savedStateId) { ++tableVersion_; });

        this->value_.fill(-std::numeric_limits<Likelihood>::infinity());
        this->setDirty();
        this->MultinomialTransmissionProcess::value();
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)> MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::value() {
        /*
         * Value is a matrix that gives the probability of transmitting some number of strains from some number of parents. The
         * rows are the number of parents, and the columns are the number of strains. We assume the number of strains transmitted
         * from a single parent is zt-poisson with rate equal to `meanStrainsTransmitted_`. In the case of multiple parents, the number
         * of strains transmitted is the sum of the number of strains transmitted from each parent
         * todo: meanStrainsTransmitted_ isn't actually the mean of the distribution
         */
        if (this->isDirty()) {
            const double lambda = meanStrainsTransmitted_->value();
            for (unsigned int kk = 0; kk < MAX_PARENTS + 1; ++kk) {
                const int num_parents = static_cast<int>(kk) + 1;
                Probability denominator = 0.0;
                for (unsigned int jj = kk; jj < MAX_STRAINS; ++jj) {
                    const int num_strains = static_cast<int>(jj) + 1;

                    double correction = 0;
                    for (int ii = 0; ii < num_parents; ++ii) {
                        correction += std::pow(-1, ii) * std::pow(num_parents - ii, num_strains) * boost::math::binomial_coefficient<double>(num_parents, ii);
                    }
                    this->value_[kk * MAX_STRAINS + jj] = num_strains * std::log(lambda) - (num_parents * std::log(std::exp(lambda) - 1)) - std::log(boost::math::factorial<double>(num_strains)) + std::log(correction);
                    denominator += std::exp(this->value_[kk * MAX_STRAINS + jj]);
                }

                // Normalize
                denominator = std::log(denominator);
                for (unsigned int jj = 0; jj < MAX_STRAINS; ++jj) {
                    this->value_[kk * MAX_STRAINS + jj] -= denominator;
                }
            }
            ++tableVersion_;
            this->setClean();
        }
        return this->value_;
    }



    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    const typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::LocusTerm&
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::locusTerm(const LocusTermKey& key) {
        if (const auto it = locusTerms_.find(key); it != locusTerms_.end()) {
            return it->second;
        }
        return locusTerms_.emplace(key, computeLocusTerm(key)).first->second;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    void MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::reserveCacheCapacity() {
        if (locusTerms_.size() >= MAX_LOCUS_TERM_CACHE_SIZE) {
            locusTerms_.clear();
        }
        if (strainSums_.size() >= MAX_STRAIN_SUM_CACHE_SIZE or genomeIds_.size() >= MAX_INTERNED_GENOMES) {
            // Strain sums are keyed by genome ids, so they are dropped with the ids.
            strainSums_.clear();
            genomeIds_.clear();
            infectionGenomes_.clear();
#ifndef NDEBUG
            genomesById_.clear();
#endif
        }
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::GenomeId
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::genomeId(const p_Infection& infection) {
        // Latent genotypes in the latentGenotype() map's order, which is the same for every infection in a chain.
        genomeScratch_.clear();
        const auto version = infection->latentGenomeVersion();
        if (const auto it = infectionGenomes_.find(infection->uid()); it != infectionGenomes_.end() and it->second.version == version) {
#ifndef NDEBUG
            for (const auto& [locus, genotype] : infection->latentGenotype()) {
                genomeScratch_.push_back(genotype->value());
            }
            assert(*genomesById_[it->second.id] == genomeScratch_ && "latent genotype changed without changing latentGenomeVersion()");
#endif
            return it->second.id;
        }

        for (const auto& [locus, genotype] : infection->latentGenotype()) {
            genomeScratch_.push_back(genotype->value());
        }
        const auto [interned, inserted] = genomeIds_.try_emplace(genomeScratch_, static_cast<GenomeId>(genomeIds_.size()));
#ifndef NDEBUG
        if (inserted) {
            genomesById_.push_back(&interned->first);
        }
#endif
        infectionGenomes_[infection->uid()] = InfectionGenome{.version = version, .id = interned->second};
        return interned->second;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::GenomeId
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::locusOrderId(const p_Infection& infection) {
        // An infection's loci are fixed once it is constructed.
        if (const auto it = infectionLocusOrders_.find(infection->uid()); it != infectionLocusOrders_.end()) {
            return it->second;
        }
        std::vector<const void*> order;
        for (const auto& locus : infection->loci()) {
            order.push_back(locus.get());
        }
        const auto [interned, inserted] = locusOrderIds_.try_emplace(std::move(order), static_cast<GenomeId>(locusOrderIds_.size()));
        infectionLocusOrders_[infection->uid()] = interned->second;
        return interned->second;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::LocusTerm
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::computeLocusTerm(const LocusTermKey& key) {
        // The arithmetic below, including operation order and the zero-probability cutoffs, matches the
        // per-variant code this replaced so that results are bit-identical.
        LocusTerm out{};
        const auto& childGenotype = key.genotypes[0];
        const size_t numParents = key.numParents;

        // Only the child's alleles' frequencies are used, so only they are accumulated, in index order; each receives
        // the same additions, in the same parent order, as when every allele was accumulated.
        childAlleles_.clear();
        for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
            if (childGenotype.allele(i)) {
                childAlleles_.push_back(i);
            }
        }
        parentPopFreqs_.assign(childAlleles_.size(), 0.0);
        double zeroProbCutoff = 1e-6;

        switch (key.variant) {
            case Variant::ObservedParents: {
                zeroProbCutoff = 1e-10;
                for (size_t p = 1; p <= numParents; ++p) {
                    const auto& parentGenotype = key.genotypes[p];
                    const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                    for (size_t k = 0; k < childAlleles_.size(); ++k) {
                        if (parentGenotype.allele(childAlleles_[k])) {
                            parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / (totalAllelesPresent * numParents));
                        }
                    }
                }
                break;
            }
            case Variant::ObservedAndLatentParents: {
                for (size_t p = 1; p < numParents; ++p) {
                    const auto& parentGenotype = key.genotypes[p];
                    const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                    for (size_t k = 0; k < childAlleles_.size(); ++k) {
                        if (parentGenotype.allele(childAlleles_[k])) {
                            parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / totalAllelesPresent / static_cast<Probability>(numParents));
                        }
                    }
                }

                const auto& latentParentGenotype = key.genotypes[numParents];
                const int totalAllelesPresent = latentParentGenotype.totalPositiveCount();
                if (GeneticsImpl::truePositiveCount(latentParentGenotype, childGenotype) == 0) {
                    out.impossible = true;
                    return out;
                }
                for (size_t k = 0; k < childAlleles_.size(); ++k) {
                    parentPopFreqs_[k] += (static_cast<Probability>(latentParentGenotype.allele(childAlleles_[k])) / totalAllelesPresent / static_cast<Probability>(numParents));
                }
                break;
            }
            case Variant::LatentParentOnly: {
                const auto& parentGenotype = key.genotypes[1];
                const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                for (size_t k = 0; k < childAlleles_.size(); ++k) {
                    if (parentGenotype.allele(childAlleles_[k])) {
                        parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / totalAllelesPresent);
                    }
                }
                break;
            }
        }

        Probability constrainedSetProb = 0.0;
        prVec_.clear();
        bool zeroProbEvent = false;

        for (size_t k = 0; k < childAlleles_.size(); ++k) {
            prVec_.push_back(parentPopFreqs_[k]);
            constrainedSetProb += parentPopFreqs_[k];
            zeroProbEvent = zeroProbEvent || std::abs(parentPopFreqs_[k]) < zeroProbCutoff;
        }

        if (prVec_.empty() || zeroProbEvent) {
            out.impossible = true;
            return out;
        }

        for (Likelihood& af : prVec_) {
            af /= constrainedSetProb;
        }

        const Likelihood logConstrainedSetProb = std::log(constrainedSetProb);
        const std::vector<Likelihood>& pamVec = probAnyMissing_.vectorized(prVec_, MAX_STRAINS);
        for (unsigned int numStrains = 1; numStrains <= MAX_STRAINS; ++numStrains) {
            const unsigned int idx = numStrains - 1;
            out.terms[idx] = pamVec[idx] >= 1.0 ? -std::numeric_limits<Likelihood>::infinity() : std::log(1.0 - pamVec[idx]) + logConstrainedSetProb * numStrains;
        }
        return out;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::StrainSumEntry&
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::strainSums(
            const Variant variant,
            const p_Infection& infection,
            const ParentSet* parentSet,
            const p_Infection& latentParent,
            const size_t numParents,
            const ResolvedGenomes* resolved) {
        StrainSumKey sumKey{.variant = variant, .numParents = static_cast<unsigned int>(numParents)};
        size_t slot = 1;
        if (resolved != nullptr) {
            sumKey.childLocusOrder = resolved->childLocusOrder;
            sumKey.genomes[0]      = resolved->child;
            if (parentSet != nullptr) {
                for (size_t i = 0; i < parentSet->size(); ++i) {
                    sumKey.genomes[slot++] = resolved->parents[i];
                }
            }
            if (latentParent != nullptr) {
                sumKey.genomes[slot] = resolved->latentParent;
            }
        } else {
            sumKey.childLocusOrder = locusOrderId(infection);
            sumKey.genomes[0]      = genomeId(infection);
            if (parentSet != nullptr) {
                for (const auto& parent : *parentSet) {
                    sumKey.genomes[slot++] = genomeId(parent);
                }
            }
            if (latentParent != nullptr) {
                sumKey.genomes[slot] = genomeId(latentParent);
            }
        }
        if (const auto it = strainSums_.find(sumKey); it != strainSums_.end()) {
            return it->second;
        }

        StrainSums sums{};
        LocusTermKey key{.variant = variant, .numParents = static_cast<unsigned int>(numParents)};
        for (const auto& locus : infection->loci()) {
            key.genotypes[0] = infection->latentGenotype(locus)->value();
            slot             = 1;
            if (parentSet != nullptr) {
                for (const auto& parent : *parentSet) {
                    key.genotypes[slot++] = parent->latentGenotype(locus)->value();
                }
            }
            if (latentParent != nullptr) {
                key.genotypes[slot] = latentParent->latentGenotype(locus)->value();
            }

            const LocusTerm& term = locusTerm(key);
            sums.add(term);
            if (term.impossible) {
                // Most parent sets of a child are impossible, mostly at many loci, so the rest are not summed.
                sums.stoppedAt = locus.get();
                break;
            }
        }
        return strainSums_.emplace(sumKey, StrainSumEntry{.sums = sums}).first->second;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    std::optional<Likelihood> MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::strainSummedLogLikelihood(
            const Variant variant,
            const p_Infection& infection,
            const ParentSet* parentSet,
            const p_Infection& latentParent,
            const size_t numParents,
            const ResolvedGenomes* resolved) {
        StrainSumEntry& entry = strainSums(variant, infection, parentSet, latentParent, numParents, resolved);
        if (entry.sums.impossible) {
            return std::nullopt;
        }
        if (this->isDirty()) {
            this->value();
        }
        if (entry.tableVersion != tableVersion_) {
            entry.summed       = strainSummedLogLikelihood(entry.sums, numParents);
            entry.tableVersion = tableVersion_;
        }
        return entry.summed;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::strainSummedLogLikelihood(const StrainSums& sums, const size_t numParents) {
        if (sums.impossible) {
            return -std::numeric_limits<Likelihood>::infinity();
        }
        if (this->isDirty()) {
            this->value();
        }
        std::array<Likelihood, MAX_STRAINS> logLikelihoods{};
        for (unsigned int idx = 0; idx < MAX_STRAINS; ++idx) {
            logLikelihoods[idx] = sums.infinite[idx] > 0 ? -std::numeric_limits<Likelihood>::infinity() : sums.finite[idx];
        }

        // Every parent transmits at least one strain, so fewer strains than parents is impossible. These entries
        // used to enter the log-sum-exp without a strain-count probability, inflating multi-parent sets.
        for (unsigned int numStrains = 1; numStrains < numParents; ++numStrains) {
            logLikelihoods[numStrains - 1] = -std::numeric_limits<Likelihood>::infinity();
        }

        Likelihood maxLl = -std::numeric_limits<Likelihood>::infinity();
        for (unsigned int numStrains = numParents; numStrains <= MAX_STRAINS; ++numStrains) {
            unsigned int idx = numStrains - 1;
            // Add the probability of the number of strains
            if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                continue;
            }
            logLikelihoods[idx] += this->probNumStrains(numStrains, numParents);
            maxLl = std::max(maxLl, logLikelihoods[idx]);
        }

        return core::utils::logSumExpKnownMax(logLikelihoods.begin(), logLikelihoods.end(), maxLl);
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    const typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::StrainSums&
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::parentSetSums(
            const Variant variant,
            const p_Infection& infection,
            const ParentSet* parentSet,
            const p_Infection& latentParent,
            const ResolvedGenomes& resolved) {
        const size_t numParents = (parentSet != nullptr ? parentSet->size() : 0) + (latentParent != nullptr ? 1 : 0);
        return strainSums(variant, infection, parentSet, latentParent, numParents, &resolved).sums;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    const typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::LocusTerm&
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::locusTerm(
            const Variant variant,
            const GeneticsImpl& child,
            const std::span<const GeneticsImpl* const> parents,
            const GeneticsImpl* latentParent) {
        LocusTermKey key{.variant = variant, .numParents = static_cast<unsigned int>(parents.size() + (latentParent != nullptr ? 1 : 0))};
        key.genotypes[0] = child;
        size_t slot      = 1;
        for (const auto* parent : parents) {
            key.genotypes[slot++] = *parent;
        }
        if (latentParent != nullptr) {
            key.genotypes[slot] = *latentParent;
        }
        return locusTerm(key);
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    typename MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::LocusTerm
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::uncachedLocusTerm(
            const Variant variant,
            const GeneticsImpl& child,
            const std::span<const GeneticsImpl* const> parents,
            const GeneticsImpl* latentParent) {
        LocusTermKey key{.variant = variant, .numParents = static_cast<unsigned int>(parents.size() + (latentParent != nullptr ? 1 : 0))};
        key.genotypes[0] = child;
        size_t slot      = 1;
        for (const auto* parent : parents) {
            key.genotypes[slot++] = *parent;
        }
        if (latentParent != nullptr) {
            key.genotypes[slot] = *latentParent;
        }
        return computeLocusTerm(key);
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::logLikelihood(
            const Variant variant,
            const p_Infection& infection,
            const ParentSet* parentSet,
            const p_Infection& latentParent,
            const p_ParentSetSizePrior& psp,
            const ResolvedGenomes* resolved) {
        // One for each observed parent, plus the latent parent if present.
        const size_t numParents = (parentSet != nullptr ? parentSet->size() : 0) + (latentParent != nullptr ? 1 : 0);
        const auto summed       = strainSummedLogLikelihood(variant, infection, parentSet, latentParent, numParents, resolved);
        if (!summed) {
            return -std::numeric_limits<Likelihood>::infinity();
        }
        Likelihood llik = *summed;

        // Add the prior on the number of parents. psp holds probabilities, not log-probabilities.
        llik += std::log(psp->value()(numParents));
        return llik;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename ParentSetSizePriorImpl, typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizePriorImpl, GeneticsImpl>::calculateLogLikelihood(p_Infection infection, const ParentSet& parentSet, p_ParentSetSizePrior psp) {
        reserveCacheCapacity();
        return logLikelihood(Variant::ObservedParents, infection, &parentSet, nullptr, psp, nullptr);
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H
