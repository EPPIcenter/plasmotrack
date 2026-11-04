//
// The Transmission memo: one per chain, the memory of Transmission kernel results keyed by genome and frequency
// identity, shared by every Infection's Transmission term. In a clonal cluster most Parent sets are made of the same
// few genomes, so locus terms, Parent-set sums and latent locus tables are kept by the identity of the genotypes and
// frequencies involved rather than recomputed per Infection.
//
// Identities: genomes are interned from their content (an Infection's id is recomputed only when its
// latentGenomeVersion() changes), each Infection's locus order likewise, and a Locus's Allele frequencies get an id so
// that the table at the frequencies a term last summed with stays available to swap out.
//
// Lifetime: a Transmission term calls beginBatch() at the start of each evaluation. It bounds every cache; ids and
// references handed out before it may be invalid after it, and generation() changes when the frequency ids are
// forgotten. Within a batch the genotypes involved must not change.
//
// See doc/marginal-latent-parent.md for the memoised layers of the latent locus table.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSMISSIONMEMO_H
#define TRANSMISSION_NETWORKS_APP_TRANSMISSIONMEMO_H

#include "core/datatypes/Genotype.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"
#include "model/transmission_process/TransmissionKernel.h"
#include "model/transmission_process/source_transmission_process/SourceOverlap.h"

#include <array>
#include <cassert>
#include <cstdint>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

namespace transmission_nets::model::transmission_process {

    template<int MAX_PARENTS, int MAX_STRAINS, int MAX_COI, typename InfectionImpl, typename LocusImpl, typename AlleleFrequencyContainer>
    class TransmissionMemo {
    public:
        using Genotype = core::datatypes::Genotype;

    public:
        using Kernel           = TransmissionKernel<MAX_PARENTS, MAX_STRAINS, MAX_COI>;
        using Likelihood       = typename Kernel::Likelihood;
        using Variant          = typename Kernel::Variant;
        using LocusTerm        = typename Kernel::LocusTerm;
        using StrainSums       = typename Kernel::StrainSums;
        using LocusTable       = typename Kernel::LocusTable;
        using LatentSums       = typename Kernel::LatentSums;
        using StrainCountTable = typename Kernel::StrainCountTable;
        static constexpr int NUM_CELLS = Kernel::NUM_CELLS;

        using GenomeId          = std::uint32_t;
        using FrequencyId       = std::uint32_t;
        using p_Infection       = std::shared_ptr<InfectionImpl>;
        using p_Locus           = std::shared_ptr<LocusImpl>;
        using ParentSet         = core::containers::ParentSet<InfectionImpl>;
        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;

        // Genome ids resolved once per evaluation rather than per Parent set.
        struct ResolvedGenomes {
            GenomeId childLocusOrder = 0;
            GenomeId child           = 0;
            GenomeId latentParent    = 0;
            std::array<GenomeId, MAX_PARENTS> parents{};// in the parent set's iteration order
        };

        TransmissionMemo(p_ParameterDouble meanStrainsTransmitted, std::shared_ptr<AlleleFrequencyContainer> alleleFrequencies)
            : meanStrainsTransmitted_(std::move(meanStrainsTransmitted)), alleleFrequencies_(std::move(alleleFrequencies)) {}

        [[nodiscard]] const p_ParameterDouble& meanStrainsTransmitted() const { return meanStrainsTransmitted_; }
        [[nodiscard]] const std::shared_ptr<AlleleFrequencyContainer>& alleleFrequencies() const { return alleleFrequencies_; }

        // Bounds every cache, each at its own limit, and only here: ids and results handed out earlier in the same
        // batch stay valid until the next call.
        void beginBatch() {
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
            if (frequencies_.size() >= MAX_FREQUENCY_IDS) {
                frequencyIds_.clear();
                frequencies_.clear();
                frequencyStore_.clear();
                lastFrequencies_.clear();
                locusStores_.clear();
                stores_.clear();
                ++generation_;
            }
        }

        // Changes when frequency ids are forgotten; ids from an earlier generation mean nothing.
        [[nodiscard]] std::uint64_t generation() const { return generation_; }

        // The id of an Infection's latent genome: Infections with the same latent genotypes share one.
        GenomeId genomeOf(const p_Infection& infection) {
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

        // The id of an Infection's locus order, which is fixed once it is constructed.
        GenomeId locusOrderOf(const p_Infection& infection) {
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

        // The id of a Locus's current Allele frequencies.
        FrequencyId frequencyId(const p_Locus& locus) {
            const auto& freqs = alleleFrequencies_->alleleFrequencies(locus)->value();
            auto& last        = lastFrequencies_[locus.get()];
            bool same         = last.valid and last.values.size() == freqs.totalElements();
            for (size_t j = 0; same and j < last.values.size(); ++j) {
                same = last.values[j] == freqs.frequencies(j);
            }
            if (same) {
                return last.id;
            }
            last.values.resize(freqs.totalElements());
            for (size_t j = 0; j < last.values.size(); ++j) {
                last.values[j] = freqs.frequencies(j);
            }
            FrequencyKey key{.locus = locus.get(), .values = last.values};
            auto [it, inserted] = frequencyIds_.try_emplace(std::move(key), static_cast<FrequencyId>(frequencies_.size()));
            if (inserted) {
                frequencies_.push_back(&it->first);
                auto [store, added] = locusStores_.try_emplace(locus.get(), static_cast<std::uint32_t>(stores_.size()));
                if (added) {
                    stores_.emplace_back();
                }
                frequencyStore_.push_back(store->second);
            }
            last.id    = it->second;
            last.valid = true;
            return last.id;
        }

        // log P(m | n), current for the strain-count parameter.
        const StrainCountTable& strainCountTable() {
            const double lambda = meanStrainsTransmitted_->value();
            if (not tableValid_ or lambda != tableLambda_) {
                table_       = Kernel::strainCountTable(lambda);
                tableLambda_ = lambda;
                tableValid_  = true;
            }
            return table_;
        }

        // One locus' term from the genotypes at it: the child's, the observed parents' in parent-set order, and the
        // latent parent's (null for the observed-parents variant). Valid until the next beginBatch().
        const LocusTerm& locusTerm(const Variant variant, const core::datatypes::Genotype& child, const std::span<const core::datatypes::Genotype* const> parents, const core::datatypes::Genotype* latentParent) {
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

        // A Parent set's locus terms summed over the child's loci. parentSet is null for the latent-parent-only variant
        // and latentParent null for the observed-parents variant. Valid until the next beginBatch().
        const StrainSums& parentSetSums(const Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, const ResolvedGenomes& resolved) {
            const size_t numParents = (parentSet != nullptr ? parentSet->size() : 0) + (latentParent != nullptr ? 1 : 0);
            return strainSums(variant, infection, parentSet, latentParent, numParents, &resolved);
        }

        // A Parent set's log-likelihood marginalised over strain counts, before the parent-set-size term. numParents
        // counts the latent parent.
        Likelihood strainSummedLogLikelihood(const StrainSums& sums, const size_t numParents) {
            return Kernel::strainSummedLogLikelihood(sums, numParents, strainCountTable());
        }

        // log g(k, m) at a locus for the child genotype and the observed parents' genotypes (none for the latent
        // parent alone), at the frequencies with the given id. Valid until the next call.
        const LocusTable& locusTable(const FrequencyId frequencyId, const core::datatypes::Genotype& child, const std::span<const core::datatypes::Genotype* const> parents) {
            auto& slot = slotFor(frequencyId);
            auto key   = genotypeKey(child, parents);
            if (const auto it = slot.tables.find(key); it != slot.tables.end()) {
                return it->second;
            }
            // Two generations: a table from the older one is moved to the current one when used, and the older one is
            // dropped when the current one fills, so tables in use survive.
            if (const auto it = slot.olderTables.find(key); it != slot.olderTables.end()) {
                auto node = slot.olderTables.extract(it);
                return slot.tables.insert(std::move(node)).position->second;
            }
            if (slot.tables.size() >= MAX_TABLES_PER_SLOT) {
                slot.olderTables = std::move(slot.tables);
                slot.tables      = {};
            }
            const auto& f = fTable(key, child, parents);
            if (f.entries.empty()) {
                auto& table      = slot.tables[std::move(key)];
                table.impossible = true;
                return table;
            }
            const marginal::SourceOverlap* overlap = nullptr;
            if (const auto it = slot.overlaps.find(child); it != slot.overlaps.end()) {
                overlap = &it->second;
            } else {
                const auto childAlleles = Kernel::allelesOf(child);
                overlap = &slot.overlaps.emplace(child, marginal::SourceOverlap(frequencies_[frequencyId]->values, childAlleles, MAX_COI)).first->second;
            }
            auto& table = slot.tables[std::move(key)];
            Kernel::contract(f, *overlap, table);
            return table;
        }

        // The Parent set's log-likelihood with the latent parent summed out, before the parent-set-size term.
        // numParents counts the latent parent.
        Likelihood latentSummedLogLikelihood(const LatentSums& sums, const size_t numParents, const std::array<double, MAX_COI + 1>& coiLogProbs) {
            if (sums.impossible) {
                return Kernel::NEG_INF;
            }
            return Kernel::latentSummedLogLikelihood(sums, numParents, coiLogProbs, strainCountTable());
        }

        // The plain term for the observed parents alone, before the parent-set-size term: the form the consistency
        // test checks the Transmission terms against.
        Likelihood observedParentsLogLikelihood(const p_Infection& infection, const ParentSet& parentSet) {
            beginBatch();
            return strainSummedLogLikelihood(strainSums(Variant::ObservedParents, infection, &parentSet, nullptr, parentSet.size(), nullptr), parentSet.size());
        }

    private:
        // A locus term depends only on these genotypes, so identical genotypes (common in clonal clusters) share one
        // evaluation. genotypes[0] is the child, then the parents in parent-set order, latent parent last.
        struct LocusTermKey {
            Variant variant;
            unsigned int numParents;
            std::array<core::datatypes::Genotype, MAX_PARENTS + 2> genotypes{};
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

        const LocusTerm& locusTerm(const LocusTermKey& key) {
            if (const auto it = locusTerms_.find(key); it != locusTerms_.end()) {
                return it->second;
            }
            // The child, the observed parents and, where the variant has one, the latent parent last.
            const unsigned int numObserved   = key.variant == Variant::ObservedParents ? key.numParents : key.numParents - 1;
            const core::datatypes::Genotype* latentParent = key.variant == Variant::ObservedParents ? nullptr : &key.genotypes[key.numParents];
            for (unsigned int p = 0; p < numObserved; ++p) {
                parentPointers_[p] = &key.genotypes[p + 1];
            }
            return locusTerms_.emplace(key, kernel_.locusTerm(key.variant, key.genotypes[0], std::span<const core::datatypes::Genotype* const>(parentPointers_.data(), numObserved), latentParent)).first->second;
        }

        // Genome-level memo. Re-evaluating a parent set means summing locus terms over every locus, but in a clonal
        // cluster most parent sets are made of the same few genomes, so the per-locus sums are cached by the identity
        // of the whole genomes involved.
        struct InfectionGenome {
            std::uint64_t version;
            GenomeId id;
        };
        struct GenomeHash {
            std::size_t operator()(const std::vector<core::datatypes::Genotype>& genome) const noexcept {
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

        // `resolved` may be null, in which case genome ids are looked up here.
        const StrainSums& strainSums(const Variant variant, const p_Infection& infection, const ParentSet* parentSet, const p_Infection& latentParent, const size_t numParents, const ResolvedGenomes* resolved) {
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
                sumKey.childLocusOrder = locusOrderOf(infection);
                sumKey.genomes[0]      = genomeOf(infection);
                if (parentSet != nullptr) {
                    for (const auto& parent : *parentSet) {
                        sumKey.genomes[slot++] = genomeOf(parent);
                    }
                }
                if (latentParent != nullptr) {
                    sumKey.genomes[slot] = genomeOf(latentParent);
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
            return strainSums_.emplace(sumKey, sums).first->second;
        }

        // The latent locus tables: see the memoised layers in doc/marginal-latent-parent.md.
        struct GenotypeKey {
            unsigned int numParents = 0;// observed parents
            std::array<core::datatypes::Genotype, MAX_PARENTS + 1> genotypes{};// child, then the observed parents
            bool operator==(const GenotypeKey&) const = default;
        };
        struct GenotypeKeyHash {
            std::size_t operator()(const GenotypeKey& key) const noexcept {
                std::size_t h = key.numParents * 0x9e3779b97f4a7c15ULL;
                for (unsigned int i = 0; i <= key.numParents; ++i) {
                    h = (h * 0x100000001b3ULL) ^ key.genotypes[i].hash();
                }
                return h;
            }
        };
        struct GeneticsHash {
            std::size_t operator()(const core::datatypes::Genotype& g) const noexcept { return g.hash(); }
        };
        // The tables and Source overlaps at one locus for one frequency id. Each locus keeps the few most recently used
        // ids: an allele-frequency proposal needs the tables at the current and the proposed frequencies, and tables at
        // frequencies no longer current are never asked for again. Tables are pure functions of their keys, so dropping
        // one only costs recomputing it.
        struct FrequencySlot {
            bool used             = false;
            FrequencyId id        = 0;
            std::uint64_t lastUse = 0;
            std::unordered_map<GenotypeKey, LocusTable, GenotypeKeyHash> tables{};
            std::unordered_map<GenotypeKey, LocusTable, GenotypeKeyHash> olderTables{};
            std::unordered_map<core::datatypes::Genotype, marginal::SourceOverlap, GeneticsHash> overlaps{};
        };
        static constexpr int SLOTS_PER_LOCUS = 4;
        using LocusStore                     = std::array<FrequencySlot, SLOTS_PER_LOCUS>;

        FrequencySlot& slotFor(const FrequencyId id) {
            auto& store = stores_[frequencyStore_[id]];
            ++useCounter_;
            FrequencySlot* oldest = &store[0];
            for (auto& slot : store) {
                if (slot.used and slot.id == id) {
                    slot.lastUse = useCounter_;
                    return slot;
                }
                if (not slot.used or (oldest->used and slot.lastUse < oldest->lastUse)) {
                    oldest = &slot;
                }
            }
            oldest->used    = true;
            oldest->id      = id;
            oldest->lastUse = useCounter_;
            oldest->tables.clear();
            oldest->olderTables.clear();
            oldest->overlaps.clear();
            return *oldest;
        }

        struct FrequencyKey {
            const LocusImpl* locus = nullptr;
            std::vector<double> values{};
            bool operator==(const FrequencyKey&) const = default;
        };
        struct FrequencyKeyHash {
            std::size_t operator()(const FrequencyKey& key) const noexcept {
                std::size_t h = std::hash<const void*>{}(key.locus);
                for (const double v : key.values) {
                    h = (h * 0x100000001b3ULL) ^ std::hash<double>{}(v);
                }
                return h;
            }
        };
        struct LastFrequencies {
            bool valid     = false;
            FrequencyId id = 0;
            std::vector<double> values{};
        };

        using FTable = typename Kernel::FTable;

        static GenotypeKey genotypeKey(const core::datatypes::Genotype& child, const std::span<const core::datatypes::Genotype* const> parents) {
            GenotypeKey key{.numParents = static_cast<unsigned int>(parents.size())};
            key.genotypes[0] = child;
            for (size_t p = 0; p < parents.size(); ++p) {
                key.genotypes[p + 1] = *parents[p];
            }
            return key;
        }

        const FTable& fTable(const GenotypeKey& key, const core::datatypes::Genotype& child, const std::span<const core::datatypes::Genotype* const> parents) {
            if (const auto it = fTables_.find(key); it != fTables_.end()) {
                return it->second;
            }
            // Two generations, as for locus tables: a table from the older one moves to the current one when used, and
            // the older one is dropped when the current one fills, so tables in use survive.
            if (const auto it = olderFTables_.find(key); it != olderFTables_.end()) {
                auto node = olderFTables_.extract(it);
                return fTables_.insert(std::move(node)).position->second;
            }
            if (fTables_.size() >= MAX_F_TABLES) {
                olderFTables_ = std::move(fTables_);
                fTables_      = {};
            }
            return fTables_.emplace(key, kernel_.fTable(child, parents)).first->second;
        }

        // Bounds on cache memory.
        static constexpr std::size_t MAX_LOCUS_TERM_CACHE_SIZE = 1 << 17;
        static constexpr std::size_t MAX_STRAIN_SUM_CACHE_SIZE = 1 << 19;
        static constexpr std::size_t MAX_INTERNED_GENOMES      = 1 << 14;
        static constexpr std::size_t MAX_TABLES_PER_SLOT       = 1 << 12;
        static constexpr std::size_t MAX_F_TABLES              = 1 << 17;
        static constexpr std::size_t MAX_FREQUENCY_IDS         = 1 << 22;

        p_ParameterDouble meanStrainsTransmitted_;
        std::shared_ptr<AlleleFrequencyContainer> alleleFrequencies_;
        Kernel kernel_{};

        // The strain-count table at the parameter's value last read.
        StrainCountTable table_{};
        double tableLambda_ = 0;
        bool tableValid_    = false;

        // Scratch.
        std::array<const core::datatypes::Genotype*, MAX_PARENTS + 1> parentPointers_{};
        std::vector<core::datatypes::Genotype> genomeScratch_{};

        // Caches are per instance (one per chain); all keys are pure functions of genotype content.
        std::unordered_map<LocusTermKey, LocusTerm, LocusTermKeyHash> locusTerms_{};
        std::unordered_map<StrainSumKey, StrainSums, StrainSumKeyHash> strainSums_{};
        std::unordered_map<std::vector<core::datatypes::Genotype>, GenomeId, GenomeHash> genomeIds_{};
        std::unordered_map<std::vector<const void*>, GenomeId, LocusOrderHash> locusOrderIds_{};
        std::unordered_map<int, InfectionGenome> infectionGenomes_{};
        std::unordered_map<int, GenomeId> infectionLocusOrders_{};
#ifndef NDEBUG
        std::vector<const std::vector<core::datatypes::Genotype>*> genomesById_{};
#endif
        std::uint64_t generation_ = 0;
        std::unordered_map<FrequencyKey, FrequencyId, FrequencyKeyHash> frequencyIds_{};
        std::vector<const FrequencyKey*> frequencies_{};// by id
        std::vector<std::uint32_t> frequencyStore_{};// by id: its locus' store
        std::unordered_map<const LocusImpl*, LastFrequencies> lastFrequencies_{};
        std::unordered_map<const LocusImpl*, std::uint32_t> locusStores_{};
        std::vector<LocusStore> stores_{};
        std::uint64_t useCounter_ = 0;
        std::unordered_map<GenotypeKey, FTable, GenotypeKeyHash> fTables_{};
        std::unordered_map<GenotypeKey, FTable, GenotypeKeyHash> olderFTables_{};
    };

}// namespace transmission_nets::model::transmission_process

#endif//TRANSMISSION_NETWORKS_APP_TRANSMISSIONMEMO_H
