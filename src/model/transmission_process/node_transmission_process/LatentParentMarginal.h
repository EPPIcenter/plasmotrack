//
// The latent parent summed out of the model. See doc/marginal-latent-parent.md.
//
// An Infection's latent parent L is drawn from the Source. A parent set with it, S ∪ {L}, enters the Infection's
// transmission term through E_{L~Source}[ P(child | S ∪ {L}) ], which factorises over loci given L's COI k and the
// number of strains transmitted m:
//   E_L[ P(child | S ∪ {L}) ] = sum_{k, m} P(k) P(m | n) prod_l g_l(k, m).
// This class gives the locus tables log g_l(k, m) and the sum over (k, m). One per chain, shared by every Infection's
// transmission process, like the node transmission process whose locus terms it builds on.
//
// Three memoised layers, each a pure function of its key:
//  - the F table: the node transmission process's locus term f(m; U, r) for every overlap U of the latent parent with
//    the child's genotype and every count r of the latent parent's alleles. Depends on the genotypes only.
//  - the Source overlap P(U, w | k) at a locus, for a child genotype. Depends on the locus' allele frequencies.
//  - the locus table log g(k, m), the two contracted.
// Allele frequencies are interned: a locus' frequencies get an id, and the frequency-dependent layers are keyed by it,
// so the table at frequencies a transmission process last summed with stays available to swap out.
//

#ifndef TRANSMISSION_NETWORKS_APP_LATENTPARENTMARGINAL_H
#define TRANSMISSION_NETWORKS_APP_LATENTPARENTMARGINAL_H

#include "core/utils/fastmath.h"
#include "model/transmission_process/TransmissionKernel.h"
#include "model/transmission_process/source_transmission_process/SourceOverlap.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

namespace transmission_nets::model::transmission_process {

    template<int MAX_PARENTS, int MAX_COI, int MAX_STRAINS, typename NodeTransmissionProcessImpl, typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GeneticsImpl, typename LocusImpl>
    class LatentParentMarginal {
    public:
        // The kernel's arithmetic and the types it sums: see TransmissionKernel.h.
        using Kernel                   = TransmissionKernel<MAX_PARENTS, MAX_STRAINS, MAX_COI, GeneticsImpl>;
        static constexpr int NUM_CELLS = Kernel::NUM_CELLS;
        using Likelihood               = double;
        using FrequencyId              = std::uint32_t;
        using Variant                  = typename NodeTransmissionProcessImpl::Variant;
        using LocusTable               = typename Kernel::LocusTable;
        using Sums                     = typename Kernel::LatentSums;

        LatentParentMarginal(std::shared_ptr<NodeTransmissionProcessImpl> ntp,
                             std::shared_ptr<COIProbabilityImpl> coiProb,
                             std::shared_ptr<AlleleFrequencyContainer> alleleFrequencies)
            : ntp_(std::move(ntp)), coiProb_(std::move(coiProb)), alleleFrequencies_(std::move(alleleFrequencies)) {}

        [[nodiscard]] const std::shared_ptr<COIProbabilityImpl>& coiProb() const { return coiProb_; }
        [[nodiscard]] const std::shared_ptr<AlleleFrequencyContainer>& alleleFrequencies() const { return alleleFrequencies_; }

        // Changes when frequency ids are forgotten; ids from an earlier generation mean nothing.
        [[nodiscard]] std::uint64_t generation() const { return generation_; }

        // The id of a locus' current allele frequencies.
        FrequencyId frequencyId(const std::shared_ptr<LocusImpl>& locus) {
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

        // log g(k, m) at a locus for the child genotype and the observed parents' genotypes (none for the latent
        // parent alone), at the frequencies with the given id. Valid until the next call.
        const LocusTable& locusTable(const FrequencyId frequencyId, const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents) {
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

        // The parent set's log-likelihood, before the parent-set-size prior, summed over the latent parent's COI and
        // the strains transmitted. numParents counts the latent parent.
        Likelihood summedLogLikelihood(const Sums& sums, const size_t numParents) {
            if (sums.impossible) {
                return Kernel::NEG_INF;
            }
            return Kernel::latentSummedLogLikelihood(sums, numParents, coiProb_->value(), ntp_->strainCountTable());
        }

        // Bounds memory. Forgets every frequency id when there are too many, which a transmission process notices
        // through generation(). Call only between evaluations, when no table reference is held.
        void reserveCapacity() {
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

    private:
        struct GenotypeKey {
            unsigned int numParents = 0;// observed parents
            std::array<GeneticsImpl, MAX_PARENTS + 1> genotypes{};// child, then the observed parents
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
            std::size_t operator()(const GeneticsImpl& g) const noexcept { return g.hash(); }
        };
        // The tables and Source overlaps at one locus for one frequency id. Each locus keeps the few most recently used
        // ids: an allele-frequency proposal needs the tables at the current and the proposed frequencies, and tables at
        // frequencies no longer current are never asked for again. Tables are pure functions of their keys, so dropping
        // one only costs recomputing it.
        struct FrequencySlot {
            bool used = false;
            FrequencyId id = 0;
            std::uint64_t lastUse = 0;
            std::unordered_map<GenotypeKey, LocusTable, GenotypeKeyHash> tables{};
            std::unordered_map<GenotypeKey, LocusTable, GenotypeKeyHash> olderTables{};
            std::unordered_map<GeneticsImpl, marginal::SourceOverlap, GeneticsHash> overlaps{};
        };
        static constexpr int SLOTS_PER_LOCUS = 4;
        using LocusStore = std::array<FrequencySlot, SLOTS_PER_LOCUS>;

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
            bool valid = false;
            FrequencyId id = 0;
            std::vector<double> values{};
        };

        using FTable = typename Kernel::FTable;

        static GenotypeKey genotypeKey(const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents) {
            GenotypeKey key{.numParents = static_cast<unsigned int>(parents.size())};
            key.genotypes[0] = child;
            for (size_t p = 0; p < parents.size(); ++p) {
                key.genotypes[p + 1] = *parents[p];
            }
            return key;
        }

        const FTable& fTable(const GenotypeKey& key, const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents) {
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

        static constexpr size_t MAX_TABLES_PER_SLOT = 1 << 12;
        static constexpr size_t MAX_F_TABLES        = 1 << 17;
        static constexpr size_t MAX_FREQUENCY_IDS   = 1 << 22;

        std::shared_ptr<NodeTransmissionProcessImpl> ntp_;
        std::shared_ptr<COIProbabilityImpl> coiProb_;
        Kernel kernel_{};
        std::shared_ptr<AlleleFrequencyContainer> alleleFrequencies_;

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

#endif//TRANSMISSION_NETWORKS_APP_LATENTPARENTMARGINAL_H
