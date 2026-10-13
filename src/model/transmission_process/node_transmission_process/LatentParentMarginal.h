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
#include "model/transmission_process/node_transmission_process/MarginalLatentParent.h"

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
        static constexpr int NUM_CELLS = (MAX_COI + 1) * MAX_STRAINS;// cell k * MAX_STRAINS + (m - 1); k = 0 is unused
        using Likelihood  = double;
        using FrequencyId = std::uint32_t;
        using Variant     = typename NodeTransmissionProcessImpl::Variant;

        // log g(k, m) at one locus. Where g(k, m) > 0 is a rectangle: k >= kMin (the fewest alleles a latent parent
        // needs for the parent set to be possible) and m >= mMin (the child's alleles). Cells outside it hold 0 and
        // are excluded by the bounds, so tables add as plain arrays. `impossible` when no cell is positive.
        struct LocusTable {
            // Leaves `terms` uninitialised: tables are built in place in the memo, which would otherwise zero them.
            LocusTable() noexcept {}
            bool impossible   = false;
            std::uint8_t kMin = 0;
            std::uint8_t mMin = 0;
            std::array<Likelihood, NUM_CELLS> terms;
        };

        // Locus tables summed over a child's loci, kept so that one locus' table can be swapped for another. The sum's
        // support is the intersection of the tables' rectangles, so the bounds are kept as counts: removing a table is
        // exact.
        struct Sums {
            std::array<Likelihood, NUM_CELLS> finite{};
            std::array<std::uint16_t, MAX_COI + 1> kMinCount{};
            std::array<std::uint16_t, MAX_STRAINS + 1> mMinCount{};
            std::uint16_t impossible = 0;
            // As StrainSums::stoppedAt: summing stopped at this locus, the first that makes the parent set impossible.
            const LocusImpl* stoppedAt = nullptr;

            void add(const LocusTable& table) {
                if (table.impossible) {
                    ++impossible;
                    return;
                }
                for (int i = 0; i < NUM_CELLS; ++i) {
                    finite[i] += table.terms[i];
                }
                ++kMinCount[table.kMin];
                ++mMinCount[table.mMin];
            }

            void remove(const LocusTable& table) {
                if (table.impossible) {
                    --impossible;
                    return;
                }
                for (int i = 0; i < NUM_CELLS; ++i) {
                    finite[i] -= table.terms[i];
                }
                --kMinCount[table.kMin];
                --mMinCount[table.mMin];
            }

            // The support's lower bounds.
            [[nodiscard]] int kLow() const {
                for (int k = MAX_COI; k > 0; --k) {
                    if (kMinCount[k] > 0) {
                        return k;
                    }
                }
                return 1;
            }
            [[nodiscard]] int mLow() const {
                for (int m = MAX_STRAINS; m > 0; --m) {
                    if (mMinCount[m] > 0) {
                        return m;
                    }
                }
                return 1;
            }
        };

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
                const auto childAlleles = allelesOf(child);
                overlap = &slot.overlaps.emplace(child, marginal::SourceOverlap(frequencies_[frequencyId]->values, childAlleles, MAX_COI)).first->second;
            }
            auto& table = slot.tables[std::move(key)];
            contract(f, *overlap, table);
            return table;
        }

        // The parent set's log-likelihood, before the parent-set-size prior, summed over the latent parent's COI and
        // the strains transmitted. numParents counts the latent parent.
        Likelihood summedLogLikelihood(const Sums& sums, const size_t numParents) {
            constexpr Likelihood negInf = -std::numeric_limits<Likelihood>::infinity();
            if (sums.impossible) {
                return negInf;
            }
            const auto coiLogProbs = coiProb_->value();
            const int mFirst       = std::max(sums.mLow(), static_cast<int>(numParents));
            const int kFirst       = sums.kLow();
            std::array<Likelihood, MAX_STRAINS + 1> strainLogProbs{};
            for (int m = mFirst; m <= MAX_STRAINS; ++m) {
                strainLogProbs[m] = ntp_->probNumStrains(m, static_cast<int>(numParents));
            }
            // The support is the rectangle k >= kFirst, m >= mFirst; the COI distribution is positive for every k >= 1.
            const int width = MAX_STRAINS - mFirst + 1;
            std::array<Likelihood, NUM_CELLS> terms;
            int count          = 0;
            Likelihood maxTerm = negInf;
            for (int k = kFirst; k <= MAX_COI; ++k) {
                const Likelihood* row = sums.finite.data() + k * MAX_STRAINS + mFirst - 1;
                const Likelihood* tx  = strainLogProbs.data() + mFirst;
                const Likelihood coi  = coiLogProbs[k];
                Likelihood* out       = terms.data() + count;
                for (int j = 0; j < width; ++j) {
                    out[j]  = coi + tx[j] + row[j];
                    maxTerm = std::max(maxTerm, out[j]);
                }
                count += width;
            }
            if (maxTerm == negInf) {
                return negInf;
            }
            // Cells more than NEGLIGIBLE below the largest add less than exp(-NEGLIGIBLE) each; usually most are.
            constexpr Likelihood NEGLIGIBLE = 40;
            std::array<Likelihood, NUM_CELLS> near;
            int numNear = 0;
            for (int i = 0; i < count; ++i) {
                const Likelihood d = terms[i] - maxTerm;
                near[numNear]      = d;
                numNear += d > -NEGLIGIBLE ? 1 : 0;
            }
            Likelihood total = 0;
            for (int i = 0; i < numNear; ++i) {
                total += core::utils::expNonPositive(near[i]);
            }
            return maxTerm + std::log(total);
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

        // f(m; U, w) for the nonempty U ⊆ C (a mask over C's alleles in index order) and w = |G \ C| up to maxW for
        // which the parent set is possible, as probabilities; only those are kept. With the bounds of where g can be
        // positive: k >= rMin, the fewest alleles of a possible latent parent, and m >= mMin.
        struct FEntry {
            std::uint32_t mask = 0;
            std::uint16_t w    = 0;
            std::uint16_t r    = 0;// |U| + w
            std::array<double, MAX_STRAINS> f{};
        };
        struct FTable {
            unsigned int maxW = 0;
            std::vector<FEntry> entries{};
            int rMin = MAX_COI + 1;
            int mMin = MAX_STRAINS + 1;
        };

        static GenotypeKey genotypeKey(const GeneticsImpl& child, std::span<const GeneticsImpl* const> parents) {
            GenotypeKey key{.numParents = static_cast<unsigned int>(parents.size())};
            key.genotypes[0] = child;
            for (size_t p = 0; p < parents.size(); ++p) {
                key.genotypes[p + 1] = *parents[p];
            }
            return key;
        }

        static std::vector<unsigned int> allelesOf(const GeneticsImpl& g) {
            std::vector<unsigned int> out;
            for (unsigned int i = 0; i < g.totalAlleles(); ++i) {
                if (g.allele(i)) {
                    out.push_back(i);
                }
            }
            return out;
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
            const auto childAlleles   = allelesOf(child);
            const unsigned int A      = child.totalAlleles();
            const auto c              = static_cast<unsigned int>(childAlleles.size());
            FTable table{};
            if (c == 0 or c > static_cast<unsigned int>(MAX_STRAINS)) {
                // No alleles, or more than the most strains transmitted can carry: impossible.
                return fTables_.emplace(key, std::move(table)).first->second;
            }
            // A child allele in no observed parent has frequency 0 unless the latent parent carries it, so masks
            // without all of them are impossible and skipped without computing their terms.
            unsigned int required = 0;
            for (unsigned int i = 0; i < c; ++i) {
                if (std::ranges::none_of(parents, [&](const GeneticsImpl* parent) { return parent->allele(childAlleles[i]); })) {
                    required |= 1u << i;
                }
            }
            table.maxW                  = std::min<unsigned int>(A - c, MAX_COI);
            const unsigned int numMasks = 1u << c;
            const Variant variant = parents.empty() ? Variant::LatentParentOnly : Variant::ObservedAndLatentParents;
            // The representative latent genotypes (marginal::representative), built incrementally: U's alleles, then
            // one more allele outside C, in index order, for each w.
            std::vector<unsigned int> outside;
            for (unsigned int a = 0; a < A and outside.size() < table.maxW; ++a) {
                if (not child.allele(a)) {
                    outside.push_back(a);
                }
            }
            GeneticsImpl latent = child;
            for (unsigned int mask = 1; mask < numMasks; ++mask) {
                if ((mask & required) != required) {
                    continue;
                }
                latent.reset();
                for (unsigned int i = 0; i < c; ++i) {
                    if (mask >> i & 1u) {
                        latent.set(childAlleles[i]);
                    }
                }
                for (unsigned int w = 0; w <= table.maxW; ++w) {
                    if (std::popcount(mask) + w > static_cast<unsigned int>(MAX_COI)) {
                        break;
                    }
                    if (w > 0) {
                        latent.set(outside[w - 1]);
                    }
                    const auto term   = ntp_->uncachedLocusTerm(variant, child, parents, &latent);
                    if (term.impossible) {
                        continue;
                    }
                    FEntry entry{.mask = mask, .w = static_cast<std::uint16_t>(w), .r = static_cast<std::uint16_t>(std::popcount(mask) + w)};
                    bool any = false;
                    for (int m = 0; m < MAX_STRAINS; ++m) {
                        entry.f[m] = std::exp(term.terms[m]);
                        any        = any or entry.f[m] > 0;
                    }
                    if (any) {
                        table.rMin = std::min(table.rMin, static_cast<int>(entry.r));
                        for (int m = 0; m < MAX_STRAINS; ++m) {
                            if (entry.f[m] > 0) {
                                table.mMin = std::min(table.mMin, m + 1);
                                break;
                            }
                        }
                        table.entries.push_back(entry);
                    }
                }
            }
            table.entries.shrink_to_fit();
            return fTables_.emplace(key, std::move(table)).first->second;
        }

        static void contract(const FTable& f, const marginal::SourceOverlap& overlap, LocusTable& out) {
            if (f.entries.empty()) {
                out.impossible = true;
                return;
            }
            // Every frequency is positive, so P(U, w | k) > 0 for k >= |U| + w, and g > 0 on k >= rMin, m >= mMin.
            const int kMin = f.rMin;
            const int mMin = f.mMin;
            std::array<double, NUM_CELLS> g;
            std::fill(g.begin() + kMin * MAX_STRAINS, g.end(), 0.0);
            for (const auto& entry : f.entries) {
                const auto& fm   = entry.f;
                const double* ps = overlap.row(entry.mask, entry.w);
                for (int k = static_cast<int>(entry.r); k <= MAX_COI; ++k) {
                    const double p = ps[k];
                    double* row    = g.data() + k * MAX_STRAINS;
                    for (int m = 0; m < MAX_STRAINS; ++m) {
                        row[m] += p * fm[m];
                    }
                }
            }
            out.impossible = false;
            out.kMin       = static_cast<std::uint8_t>(kMin);
            out.mMin       = static_cast<std::uint8_t>(mMin);
            std::fill(out.terms.begin(), out.terms.end(), 0.0);
            // Inside the rectangle every cell is positive in exact arithmetic; one that underflowed below the smallest
            // normal double counts as that, which is as good as zero.
            for (int k = kMin; k <= MAX_COI; ++k) {
                const double* in = g.data() + k * MAX_STRAINS + mMin - 1;
                double* to       = out.terms.data() + k * MAX_STRAINS + mMin - 1;
                for (int j = 0; j <= MAX_STRAINS - mMin; ++j) {
                    to[j] = core::utils::logPositive(std::max(in[j], std::numeric_limits<double>::min()));
                }
            }
        }

        static constexpr size_t MAX_TABLES_PER_SLOT = 1 << 12;
        static constexpr size_t MAX_F_TABLES        = 1 << 17;
        static constexpr size_t MAX_FREQUENCY_IDS   = 1 << 22;

        std::shared_ptr<NodeTransmissionProcessImpl> ntp_;
        std::shared_ptr<COIProbabilityImpl> coiProb_;
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
