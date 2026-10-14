//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H


#include "core/computation/OrderDerivedParentSet.h"
#include "core/computation/PartialLikelihood.h"
#include "core/config.h"
#include "core/containers/Infection.h"
#include "core/io/serialize.h"
#include "core/utils/numerics.h"


#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/range/adaptors.hpp>

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <optional>
#include <span>
#include <vector>


namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<typename InfectionEventImpl>
    struct ParentSetDist {
        Likelihood sourceLlik = 0;
        std::vector<std::pair<Likelihood, core::containers::ParentSet<InfectionEventImpl>>> parentSetLliks{};
        Likelihood totalLlik = 0;
    };

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    class OrderBasedTransmissionProcessV3 : public core::computation::PartialLikelihood {
        static_assert(ParentSetMaxCardinality == 1 or ParentSetMaxCardinality == 2, "parent sets are enumerated for at most two observed parents");

        /*
         * The order based transmission process considers the set of all possible parent sets under the given ordering.
         * Assumes each infection has a potential latent unobserved parent drawn from the Source. The latent parent is
         * summed out (the Transmission memo's latent locus tables): a parent set with it has likelihood
         * E_{L~Source}[ P(child | S ∪ {L}) ],
         * and latentParent_ only names it in parent-set distributions.
         */

        using ListenerIdMap = boost::container::flat_map<std::string, core::abstract::ListenerId_t>;

    public:
        OrderBasedTransmissionProcessV3(std::shared_ptr<MemoImpl> memo,
                                        std::shared_ptr<COIProbabilityImpl> coiProb,
                                        std::shared_ptr<ParentSetSizeLikelihoodImpl> psp,
                                        std::shared_ptr<InfectionEventImpl> child,
                                        std::shared_ptr<ParentSetImpl> parent_set,
                                        std::shared_ptr<InfectionEventImpl> latent_parent,
                                        bool null_model = false);

        Likelihood value() override;

        // Allows querying the individual parent set likelihoods
        ParentSetDist<InfectionEventImpl> calcParentSetDist();

        std::string identifier() override;

        Likelihood peek() noexcept override;

        // Input parameters
        std::shared_ptr<MemoImpl> memo_;
        std::shared_ptr<COIProbabilityImpl> coiProb_;
        std::shared_ptr<ParentSetSizeLikelihoodImpl> psp_;
        std::shared_ptr<InfectionEventImpl> child_;
        std::shared_ptr<ParentSetImpl> parentSet_;
        std::shared_ptr<InfectionEventImpl> latentParent_;

        bool null_model_ = false;

    private:
        void postSaveState(int savedStateId);

        void postAcceptState();

        void postRestoreState(int savedStateId);

        using Genetics   = typename InfectionEventImpl::GeneticsType;
        using LocusPtr   = std::shared_ptr<typename InfectionEventImpl::LocusType>;
        using StrainSums  = typename MemoImpl::StrainSums;
        using LatentSums  = typename MemoImpl::LatentSums;
        using FrequencyId = typename MemoImpl::FrequencyId;
        using Variant     = typename MemoImpl::Variant;
        using GenomeId    = typename MemoImpl::GenomeId;
        using ParentSet  = core::containers::ParentSet<InfectionEventImpl>;

        // Lowers every entry's state, or only the latent-parent entries', to at most `state`.
        void degrade(unsigned char state, bool latentOnly);
        // Called before a latent genotype of the child or a parent changes.
        void recordChange(const InfectionEventImpl* infection, const LocusPtr& locus, const Genetics& previous);
        // A latent genotype that changed since the sums were last brought up to date, with its value then.
        struct GenotypeChange {
            const InfectionEventImpl* infection;
            LocusPtr locus;
            Genetics previous;
        };

        // A parent enters the likelihood only through its latent genotypes, so parent sets whose parents have the same
        // genomes have the same likelihood. The parents are partitioned into classes of identical genomes, and one entry
        // per combination of classes stands for every parent set it covers: n_c single parents from class c, n_c n_d
        // pairs from classes c < d and n_c (n_c - 1) / 2 pairs within class c. In a clonal cluster most candidate parents
        // share one genome, so an infection with dozens of candidate parents has a handful of entries, not thousands.
        //
        // Entries are stored by position, so lookups are array indexing. Position 0 is the latent parent alone. Class
        // combination g is at 1 + 2g without the latent parent and 2 + 2g with it, combinations being numbered class by
        // class: class k alone, then its pairs with classes 0, ..., k. A class added last adds entries after all others.
        //
        // Each entry keeps its locus terms summed over the child's loci. A proposal changes one locus at a time, so
        // instead of re-summing every locus, the changed loci's old terms are swapped for new ones; `changes` records
        // what changed since the sums were brought up to date. Being part of the cache, it is saved and restored
        // with the sums it refers to. A parent whose genome changes leaves its class; if no class has its new genome,
        // it starts one whose entries begin as those of the class it left, and the changed loci are then swapped.
        //
        // The log-sum-exp is taken over two groups of entries, those without the latent parent and those with it; the
        // source term every entry of the second group shares is added to its sum rather than to each entry, so a change
        // to the source process (an allele frequency, say) recomputes no entry and copies no cache. An entry's term is
        // its log-likelihood plus the log of the number of parent sets it stands for. Each group keeps its entries'
        // weights[i] = exp(term_i - weightMax); its maximum rarely changes between evaluations, so only entries changed
        // since the last one need a new exp.
        struct LogSumExpGroup {
            Likelihood weightMax = std::numeric_limits<Likelihood>::quiet_NaN();
            size_t weightArgmax  = NO_POSITION;// an entry whose term equals weightMax, if known
            Likelihood logSum    = -std::numeric_limits<Likelihood>::infinity();
        };
        static constexpr size_t NO_POSITION = std::numeric_limits<size_t>::max();
        struct LikelihoodCache {
            // The parent set the classes partition, in ParentSet order, and each parent's class.
            std::vector<const InfectionEventImpl*> parents{};
            std::vector<std::uint32_t> parentClass{};
            // Each class's size, and a member (an index into `parents`) whose genotypes are the class's.
            std::vector<std::uint32_t> classSize{};
            std::vector<std::uint32_t> classMember{};
            std::vector<StrainSums> sums{};// entries without the latent parent (odd positions)
            std::vector<LatentSums> latentSums{};// entries with it, position p at p / 2
            // The id of each child locus' allele frequencies the latent entries were summed with, in child loci order,
            // and the generation of the ids.
            std::vector<FrequencyId> frequencyIds{};
            std::uint64_t frequencyGeneration = 0;
            std::vector<std::uint16_t> updates{};// locus swaps since the sums were last computed from scratch
            std::vector<Likelihood> summed{};// strain-summed log-likelihood, before the parent-set-size term
            std::vector<Likelihood> lliks{};// one parent set's
            std::vector<Likelihood> logMultiplicity{};// log of the number of parent sets the entry stands for
            std::vector<unsigned char> valid{};// one of the entry states below
            std::vector<Likelihood> weights{};// current for VALID entries, when their group's weightMax is not NaN
            // Entries without the latent parent (odd positions), then those with it (position 0 and even positions).
            std::array<LogSumExpGroup, 2> groups{};
            std::vector<GenotypeChange> changes{};
            // Something changed since the groups' sums were brought up to date.
            bool stale = true;
        };
        // Entry states, each needing less work than the one before.
        static constexpr unsigned char INVALID       = 0;// sums must be computed from scratch; also any entry standing for no parent set, which is not kept up to date
        static constexpr unsigned char SUMS_CHANGED  = 1;// sums current, strain-summed likelihood not
        static constexpr unsigned char OUTER_CHANGED = 2;// strain-summed likelihood current, parent-set-size term not
        static constexpr unsigned char COMPUTED      = 3;// likelihood current, weight not yet updated
        static constexpr unsigned char VALID         = 4;
        // Swapping terms in and out accumulates rounding, so sums are recomputed after this many swaps.
        static constexpr std::uint16_t MAX_LOCUS_UPDATES = 1024;
        static constexpr std::uint32_t NO_CLASS          = std::numeric_limits<std::uint32_t>::max();

        // Partitions the current parent set by genome and re-keys the cache to the new classes.
        void regroup(const ParentSet& ps, bool parentSetMayHaveChanged);
        // regroup() when only `changed`'s genome changed and the parent set did not: it joins another class, keeps its
        // own, or starts one. False if that would move entries, which regroup() then does.
        bool regroupOne(const ParentSet& ps, const InfectionEventImpl* changed);
        // Resizes every entry vector for this many classes; new entries are INVALID.
        static void resizeEntries(LikelihoodCache& cache, size_t numClasses);
        static void copyEntry(LikelihoodCache& to, size_t toPosition, const LikelihoodCache& from, size_t fromPosition);
        // Brings each entry's multiplicity up to date with the class sizes.
        static void updateMultiplicities(LikelihoodCache& cache);
        // Swaps the recorded loci's terms in every affected entry's sums for their current ones.
        void applyChanges(LikelihoodCache& current);
        // Recomputes every latent entry if its frequency ids are from an earlier generation.
        void syncFrequencyGeneration(LikelihoodCache& current);
        // Swaps the locus tables of loci whose allele frequencies changed in every latent entry's sums.
        void applyFrequencyChanges(LikelihoodCache& current);
        // A latent entry's sums over the child's loci, from scratch, for the observed parents given (none for the
        // latent parent alone).
        LatentSums latentSumsFromScratch(const LikelihoodCache& current, std::span<const InfectionEventImpl* const> parents);
        // Brings a group's weights and sum up to date; its entries are at positions first, first + 2, ...
        static void updateGroup(LikelihoodCache& current, LogSumExpGroup& group, size_t first);
        // The node's log-likelihood from the groups' sums.
        static Likelihood combine(const LikelihoodCache& current);

        static size_t numEntries(size_t numClasses);
        // g for a combination of classes, in ascending order (a pair may repeat a class).
        static size_t combinationIndex(std::span<const unsigned int> classes);

        // Calls visit(classes) for every combination of up to ParentSetMaxCardinality classes, a pair possibly within
        // one class, in position order (g = 0, 1, ...).
        template<typename Visit>
        static void forEachCombination(const size_t numClasses, Visit&& visit) {
            std::array<unsigned int, 2> indices{};
            for (unsigned int k = 0; k < numClasses; ++k) {
                indices[0] = k;
                visit(std::span<const unsigned int>(indices.data(), 1));
                if constexpr (ParentSetMaxCardinality == 2) {
                    indices[1] = k;
                    for (unsigned int c = 0; c <= k; ++c) {
                        indices[0] = c;
                        visit(std::span<const unsigned int>(indices.data(), 2));
                    }
                }
            }
        }

        // The cache as of the current saved state.
        const LikelihoodCache& cache() const {
            return caches_[slots_[cacheIndex_]];
        }
        // The cache to change. Copies it into the current saved state's own slot first if it still shares one.
        LikelihoodCache& mutableCache() {
            if (slots_[cacheIndex_] != cacheIndex_) {
                // Copy-assignment reuses the destination's capacity, so this does not allocate once warmed up.
                caches_[cacheIndex_]  = caches_[slots_[cacheIndex_]];
                slots_[cacheIndex_] = cacheIndex_;
            }
            return caches_[cacheIndex_];
        }

        // One cache per nested checkpoint. Most saved states are accepted or restored without this process's cache
        // changing, so a saved state shares the slot of the one below it until its first change; slots_[i] is the
        // slot that holds saved state i's cache, never above i.
        std::array<LikelihoodCache, core::config::MAX_CACHE_SIZE> caches_{};
        std::array<size_t, core::config::MAX_CACHE_SIZE> slots_{};
        size_t cacheIndex_ = 0;
        // Scratch for value(): each parent's genome id in the current call, resolved when first needed.
        std::vector<std::optional<GenomeId>> parentGenomes_{};
        // Scratch for regroup().
        struct NewClass {
            std::uint32_t source = NO_CLASS;// the class whose entries this one's begin as, if any
            std::uint32_t member = 0;
            std::uint32_t size   = 0;
            std::optional<GenomeId> genome{};
        };
        std::vector<const InfectionEventImpl*> changedParents_{};
        std::vector<std::uint32_t> oldIndex_{};
        std::vector<std::uint32_t> continuation_{};
        std::vector<std::uint32_t> newParentClass_{};
        std::vector<NewClass> newClasses_{};
        std::vector<std::uint32_t> classOrder_{};
        std::vector<std::uint32_t> classRank_{};
        // Scratch for applyChanges(): each class's genotype at one locus before and after, and whether it changed.
        std::vector<Genetics> previousClasses_{};
        std::vector<Genetics> currentClasses_{};
        std::vector<unsigned char> classChanged_{};
        std::vector<std::pair<const InfectionEventImpl*, std::uint32_t>> changedMembers_{};
        std::vector<LocusPtr> changedLoci_{};
        // The parent set may differ from the one the cache refers to. Not checkpointed: at worst it causes a full
        // evaluation that finds nothing to do.
        bool parentsMayHaveChanged_ = true;
        // Likewise for the allele frequencies at the child's loci, and the loci where they may have.
        bool frequenciesMayHaveChanged_ = true;
        std::vector<size_t> pendingFrequencyLoci_{};
        std::vector<unsigned char> frequencyPending_{};
        // Each child locus' position in child_->loci().
        boost::container::flat_map<const typename InfectionEventImpl::LocusType*, size_t> locusIndex_{};
        // Scratch: the observed parents' genotypes at one locus.
        std::array<const Genetics*, ParentSetMaxCardinality> parentGenotypes_{};

        std::string lastUpdated_ = "None";
    };

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::OrderBasedTransmissionProcessV3(
            std::shared_ptr<MemoImpl> memo,
            std::shared_ptr<COIProbabilityImpl> coiProb,
            std::shared_ptr<ParentSetSizeLikelihoodImpl> psp,
            std::shared_ptr<InfectionEventImpl> child,
            std::shared_ptr<ParentSetImpl> parent_set,
            std::shared_ptr<InfectionEventImpl> latent_parent,
            const bool null_model)
        : memo_(std::move(memo)),
          coiProb_(std::move(coiProb)),
          psp_(std::move(psp)),
          child_(std::move(child)),
          parentSet_(std::move(parent_set)),
          latentParent_(std::move(latent_parent)),
          null_model_(null_model) {

        // The strain-count parameter enters every entry's strain-summed likelihood, through the memo's table of it.
        memo_->meanStrainsTransmitted()->add_post_change_listener([=, this]() {
            lastUpdated_ = "strain count parameter updated";
            degrade(SUMS_CHANGED, false);
            setDirty();
        });
        memo_->meanStrainsTransmitted()->registerCacheableCheckpointTarget(this);

        // The latent parent's COI distribution enters the latent entries' sum over its COI. Like the parent-set-size
        // prior below, it is read in every evaluation.
        coiProb_->add_set_dirty_listener([=, this]() {
            lastUpdated_ = "coi distribution updated";
            degrade(SUMS_CHANGED, true);
            setDirty();
        });
        coiProb_->registerCacheableCheckpointTarget(this);

        // The allele frequencies at the child's loci enter the latent entries' locus tables. A change is found by
        // comparing frequency ids at the next evaluation.
        frequencyPending_.assign(child_->loci().size(), 0);
        for (size_t l = 0; l < child_->loci().size(); ++l) {
            const auto& locus = child_->loci()[l];
            auto frequencies  = memo_->alleleFrequencies()->alleleFrequencies(locus);
            frequencies->add_post_change_listener([=, this]() {
                lastUpdated_               = "allele frequencies updated";
                frequenciesMayHaveChanged_ = true;
                if (not frequencyPending_[l]) {
                    frequencyPending_[l] = 1;
                    pendingFrequencyLoci_.push_back(l);
                }
                setDirty();
            });
            frequencies->registerCacheableCheckpointTarget(this);
            locusIndex_[locus.get()] = l;
        }

        // The parent-set-size prior enters every entry's outer term. Its value is only recomputed when read, so
        // this process must be evaluated, reading it, before the change is accepted: an accept marks it clean.
        psp_->add_set_dirty_listener([=, this]() {
            lastUpdated_ = "psp updated";
            degrade(OUTER_CHANGED, false);
            setDirty();
        });
        psp_->registerCacheableCheckpointTarget(this);

        // Genotype changes are recorded before they happen, with the genotype replaced, and applied to the sums
        // at the next evaluation.
        child_->add_latent_genotype_pre_change_listener([=, this](const auto& locus, const auto& previous) {
            recordChange(child_.get(), locus, previous);
        });
        child_->add_post_change_listener([=, this]() {
            lastUpdated_ = "child updated";
            setDirty();
        });
        child_->registerCacheableCheckpointTarget(this);

        // Saved, accepted and restored with every parent the parent set could hold, whose genotype changes the
        // cache records.
        for (const auto& parent : parentSet_->allowedParents_) {
            parent->add_latent_genotype_pre_change_listener([=, this, parent = parent.get()](const auto& locus, const auto& previous) {
                recordChange(parent, locus, previous);
            });
            parent->registerCacheableCheckpointTarget(this);
        }

        parentSet_->add_element_changed_listener([=, this]([[maybe_unused]] const auto& parent) {
            lastUpdated_ = "parent updated";
            setDirty();
        });

        // Parents joining or leaving the parent set change the classes, which the next evaluation recomputes.
        parentSet_->add_element_added_listener([=, this]([[maybe_unused]] const auto& parent) {
            lastUpdated_           = "parent added";
            parentsMayHaveChanged_ = true;
            setDirty();
        });

        parentSet_->add_element_removed_listener([=, this]([[maybe_unused]] const auto& parent) {
            lastUpdated_           = "parent removed";
            parentsMayHaveChanged_ = true;
            setDirty();
        });

        parentSet_->registerCacheableCheckpointTarget(this);

        this->addPostSaveHook([=, this](const auto savedStateID) { this->postSaveState(savedStateID); });
        this->addPostRestoreHook([=, this](const auto savedStateID) { this->postRestoreState(savedStateID); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        this->value_ = -std::numeric_limits<Likelihood>::infinity();
        this->setDirty();
        this->OrderBasedTransmissionProcessV3::value();
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    size_t OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::numEntries(const size_t numClasses) {
        size_t combinations = numClasses;
        if constexpr (ParentSetMaxCardinality == 2) {
            combinations += numClasses * (numClasses + 1) / 2;
        }
        return 1 + 2 * combinations;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    size_t OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::combinationIndex(const std::span<const unsigned int> classes) {
        if constexpr (ParentSetMaxCardinality == 1) {
            return classes[0];
        }
        // Each class j before the last one in the combination, k, has j + 2 combinations: itself, then its pairs.
        const size_t k     = classes.back();
        const size_t block = k * (k + 3) / 2;
        return classes.size() == 1 ? block : block + 1 + classes[0];
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::resizeEntries(LikelihoodCache& cache, const size_t numClasses) {
        const size_t size = numEntries(numClasses);
        cache.sums.resize(size);
        cache.latentSums.resize((size + 1) / 2);
        cache.updates.resize(size, 0);
        cache.summed.resize(size, 0.0);
        cache.lliks.resize(size, 0.0);
        cache.logMultiplicity.resize(size, std::numeric_limits<Likelihood>::quiet_NaN());
        cache.valid.resize(size, INVALID);
        cache.weights.resize(size, 0.0);
        for (auto& group : cache.groups) {
            if (group.weightArgmax != NO_POSITION and group.weightArgmax >= size) {
                group.weightArgmax = NO_POSITION;
            }
        }
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::copyEntry(LikelihoodCache& to, const size_t toPosition, const LikelihoodCache& from, const size_t fromPosition) {
        if (toPosition % 2 == 0) {
            to.latentSums[toPosition / 2] = from.latentSums[fromPosition / 2];
        } else {
            to.sums[toPosition] = from.sums[fromPosition];
        }
        to.updates[toPosition]         = from.updates[fromPosition];
        to.summed[toPosition]          = from.summed[fromPosition];
        to.lliks[toPosition]           = from.lliks[fromPosition];
        to.logMultiplicity[toPosition] = from.logMultiplicity[fromPosition];
        to.valid[toPosition]           = from.valid[fromPosition];
        to.weights[toPosition]         = from.weights[fromPosition];
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    bool OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::regroupOne(const ParentSet& ps, const InfectionEventImpl* changed) {
        const auto& current = cache();
        const auto found    = std::lower_bound(current.parents.begin(), current.parents.end(), changed, core::containers::UidLess{});
        if (found == current.parents.end() or *found != changed) {
            return false;
        }
        const auto i            = static_cast<std::uint32_t>(found - current.parents.begin());
        const auto from         = current.parentClass[i];
        const size_t numClasses = current.classSize.size();
        const bool alone        = current.classSize[from] == 1;

        // Another class with the parent's new genome. Their members are unchanged, so have the genomes summed.
        const auto genome = memo_->genomeOf(ps.begin()[i]);
        std::uint32_t to  = NO_CLASS;
        for (std::uint32_t k = 0; k < numClasses and to == NO_CLASS; ++k) {
            if (k != from and memo_->genomeOf(ps.begin()[current.classMember[k]]) == genome) {
                to = k;
            }
        }
        if (to == NO_CLASS and alone) {
            return true;// keeps its class, whose entries applyChanges() brings up to date
        }
        if (to != NO_CLASS and alone and from + 1 != numClasses) {
            return false;// emptying a class other than the last would move entries
        }

        auto& target = mutableCache();
        // A class the parent leaves takes another member if it was the one standing in for it.
        auto keepMember = [&]() {
            if (target.classMember[from] == i) {
                for (std::uint32_t j = 0; j < target.parents.size(); ++j) {
                    if (j != i and target.parentClass[j] == from) {
                        target.classMember[from] = j;
                        break;
                    }
                }
            }
        };
        if (to != NO_CLASS) {
            target.parentClass[i] = to;
            ++target.classSize[to];
            if (alone) {
                // Its class, the last, empties.
                target.classSize.pop_back();
                target.classMember.pop_back();
                resizeEntries(target, numClasses - 1);
            } else {
                --target.classSize[from];
                keepMember();
            }
        } else {
            // Starts a class of its own, whose entries begin as those of the class it left.
            --target.classSize[from];
            keepMember();
            const auto added      = static_cast<unsigned int>(numClasses);
            target.parentClass[i] = added;
            target.classSize.push_back(1);
            target.classMember.push_back(i);
            resizeEntries(target, numClasses + 1);
            std::array<unsigned int, 2> source{from, 0};
            auto copyCombination = [&](const std::span<const unsigned int> combination) {
                const size_t to = 1 + 2 * combinationIndex(combination);
                const size_t of = 1 + 2 * combinationIndex(std::span<const unsigned int>(source.data(), combination.size()));
                copyEntry(target, to, target, of);
                copyEntry(target, to + 1, target, of + 1);
            };
            std::array<unsigned int, 2> combination{added, 0};
            copyCombination(std::span<const unsigned int>(combination.data(), 1));
            if constexpr (ParentSetMaxCardinality == 2) {
                // Its pair within itself stands for no parent set while it has one member, so is left INVALID.
                for (unsigned int c = 0; c < added; ++c) {
                    combination = {c, added};
                    source      = {std::min(c, from), std::max(c, from)};
                    copyCombination(std::span<const unsigned int>(combination.data(), 2));
                }
            }
        }
        updateMultiplicities(target);
        return true;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::regroup(const ParentSet& ps, const bool parentSetMayHaveChanged) {
        const auto& current       = cache();
        const size_t n            = ps.size();
        const size_t oldClasses   = current.classSize.size();
        auto parentAt             = [&](const size_t i) -> const std::shared_ptr<InfectionEventImpl>& { return ps.begin()[i]; };

        // Parents whose genotypes differ somewhere from those last summed. Every other parent keeps its class.
        changedParents_.clear();
        for (const auto& change : current.changes) {
            if (change.infection != child_.get() and
                change.infection->latentGenotype(change.locus)->value() != change.previous and
                std::ranges::find(changedParents_, change.infection) == changedParents_.end()) {
                changedParents_.push_back(change.infection);
            }
        }
        if (not parentSetMayHaveChanged and current.parents.size() == ps.size() and not current.valid.empty()) {
            if (changedParents_.empty()) {
                return;// every parent keeps its class
            }
            if (changedParents_.size() == 1 and regroupOne(ps, changedParents_.front())) {
                return;
            }
        }

        // Each current parent's index in the parent set the classes partition. Both are in ParentSet (uid) order.
        oldIndex_.assign(n, NO_CLASS);
        for (size_t i = 0, j = 0; i < n and j < current.parents.size();) {
            const auto* parent = parentAt(i).get();
            if (parent == current.parents[j]) {
                oldIndex_[i++] = j++;
            } else if (core::containers::UidLess{}(parent, current.parents[j])) {
                ++i;
            } else {
                ++j;
            }
        }

        // A class continues while an unchanged member is still a parent.
        newClasses_.clear();
        continuation_.assign(oldClasses, NO_CLASS);
        newParentClass_.assign(n, NO_CLASS);
        for (size_t i = 0; i < n; ++i) {
            if (oldIndex_[i] == NO_CLASS or std::ranges::find(changedParents_, parentAt(i).get()) != changedParents_.end()) {
                continue;
            }
            const auto c = current.parentClass[oldIndex_[i]];
            if (continuation_[c] == NO_CLASS) {
                continuation_[c] = newClasses_.size();
                newClasses_.push_back(NewClass{.source = c, .member = static_cast<std::uint32_t>(i)});
            }
            newParentClass_[i] = continuation_[c];
            ++newClasses_[continuation_[c]].size;
        }
        // Every other parent joins the class with its genome, or starts one. A parent that changed genome starts its
        // class from the entries of the class it left, which applyChanges() brings up to date.
        for (size_t i = 0; i < n; ++i) {
            if (newParentClass_[i] != NO_CLASS) {
                continue;
            }
            const auto genome    = memo_->genomeOf(parentAt(i));
            std::uint32_t joined = NO_CLASS;
            for (std::uint32_t k = 0; k < newClasses_.size() and joined == NO_CLASS; ++k) {
                auto& candidate = newClasses_[k];
                if (not candidate.genome) {
                    candidate.genome = memo_->genomeOf(parentAt(candidate.member));
                }
                if (*candidate.genome == genome) {
                    joined = k;
                }
            }
            if (joined == NO_CLASS) {
                joined = newClasses_.size();
                newClasses_.push_back(NewClass{.source = oldIndex_[i] == NO_CLASS ? NO_CLASS : current.parentClass[oldIndex_[i]],
                                               .member = static_cast<std::uint32_t>(i),
                                               .genome = genome});
            }
            newParentClass_[i] = joined;
            ++newClasses_[joined].size;
        }

        // Ordered by the class each begins as, so that when every class has exactly one successor, as when a parent
        // with a genome of its own changes it, the entries stay where they are.
        const size_t numClasses = newClasses_.size();
        classOrder_.resize(numClasses);
        std::iota(classOrder_.begin(), classOrder_.end(), 0);
        std::ranges::stable_sort(classOrder_, std::less{}, [&](const std::uint32_t k) { return newClasses_[k].source; });
        classRank_.resize(numClasses);
        for (size_t k = 0; k < numClasses; ++k) {
            classRank_[classOrder_[k]] = k;
        }
        bool inPlace = not current.valid.empty() and numClasses == oldClasses;
        for (size_t k = 0; inPlace and k < numClasses; ++k) {
            inPlace = newClasses_[classOrder_[k]].source == k;
        }

        auto fill = [&](LikelihoodCache& target) {
            target.parents.resize(n);
            target.parentClass.resize(n);
            for (size_t i = 0; i < n; ++i) {
                target.parents[i]     = parentAt(i).get();
                target.parentClass[i] = classRank_[newParentClass_[i]];
            }
            target.classSize.resize(numClasses);
            target.classMember.resize(numClasses);
            for (size_t k = 0; k < numClasses; ++k) {
                target.classSize[k]   = newClasses_[classOrder_[k]].size;
                target.classMember[k] = newClasses_[classOrder_[k]].member;
            }
        };

        if (inPlace) {
            auto& target = mutableCache();
            fill(target);
            updateMultiplicities(target);
            return;
        }

        const size_t expected = numEntries(numClasses);
        LikelihoodCache next;
        fill(next);
        next.sums.assign(expected, StrainSums{});
        next.latentSums.assign((expected + 1) / 2, LatentSums{});
        next.frequencyIds        = current.frequencyIds;
        next.frequencyGeneration = current.frequencyGeneration;
        next.updates.assign(expected, 0);
        next.summed.assign(expected, 0.0);
        next.lliks.assign(expected, 0.0);
        next.logMultiplicity.assign(expected, std::numeric_limits<Likelihood>::quiet_NaN());
        next.valid.assign(expected, INVALID);
        next.weights.assign(expected, 0.0);
        next.groups = current.groups;
        for (auto& group : next.groups) {
            group.weightArgmax = NO_POSITION;
        }
        next.changes = current.changes;
        next.stale   = true;
        auto carry   = [&](const size_t to, const size_t from) { copyEntry(next, to, current, from); };

        if (not current.lliks.empty()) {
            // The latent-parent-only entry does not depend on the observed parents.
            carry(0, 0);

            std::array<unsigned int, 2> sources{};
            size_t position = 1;
            forEachCombination(numClasses, [&](const std::span<const unsigned int> combination) {
                bool known = true;
                for (size_t t = 0; t < combination.size(); ++t) {
                    sources[t] = newClasses_[classOrder_[combination[t]]].source;
                    known      = known and sources[t] != NO_CLASS;
                }
                if (known) {
                    // Classes are ordered by source, so a pair's sources are in order too.
                    assert(combination.size() == 1 or sources[0] <= sources[1]);
                    const size_t from = 1 + 2 * combinationIndex(std::span<const unsigned int>(sources.data(), combination.size()));
                    carry(position, from);
                    carry(position + 1, from + 1);
                }
                position += 2;
            });
        }
        updateMultiplicities(next);
        // Replaces the whole cache, so there is nothing to copy from a shared slot.
        caches_[cacheIndex_] = std::move(next);
        slots_[cacheIndex_]  = cacheIndex_;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::updateMultiplicities(LikelihoodCache& cache) {
        auto set = [&](const size_t position, const Likelihood logMultiplicity) {
            if (cache.logMultiplicity[position] == logMultiplicity) {
                return;
            }
            cache.logMultiplicity[position] = logMultiplicity;
            cache.stale                     = true;
            if (logMultiplicity == -std::numeric_limits<Likelihood>::infinity()) {
                cache.valid[position] = INVALID;// stands for no parent set, so is not kept up to date
            } else if (cache.valid[position] == VALID) {
                cache.valid[position] = COMPUTED;// its weight changes
            }
        };
        set(0, 0.0);
        size_t position = 1;
        forEachCombination(cache.classSize.size(), [&](const std::span<const unsigned int> combination) {
            const double size = cache.classSize[combination[0]];
            double count      = size;
            if (combination.size() == 2) {
                count = combination[0] == combination[1] ? size * (size - 1) / 2 : size * cache.classSize[combination[1]];
            }
            const Likelihood logMultiplicity = std::log(count);
            set(position, logMultiplicity);
            set(position + 1, logMultiplicity);
            position += 2;
        });
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::degrade(const unsigned char state, const bool latentOnly) {
        if (cache().valid.empty()) {
            return;
        }
        auto& current = mutableCache();
        current.stale = true;
        // The latent-parent entries are position 0 and every even position.
        const size_t step = latentOnly ? 2 : 1;
        current.valid[0]  = std::min(current.valid[0], state);
        for (size_t position = latentOnly ? 2 : 1; position < current.valid.size(); position += step) {
            current.valid[position] = std::min(current.valid[position], state);
        }
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::recordChange(const InfectionEventImpl* infection, const LocusPtr& locus, const Genetics& previous) {
        const auto& current = cache();
        if (null_model_ or current.valid.empty()) {
            return;
        }
        // A parent that is not in the parent set the entries refer to has no entries to update; any it gains are
        // computed from scratch.
        if (infection != child_.get() and
            not std::binary_search(current.parents.begin(), current.parents.end(), infection, core::containers::UidLess{})) {
            return;
        }
        // Only the first change since the sums were brought up to date matters: that genotype is the one summed.
        for (const auto& change : current.changes) {
            if (change.infection == infection and change.locus == locus) {
                return;
            }
        }
        auto& changed = mutableCache();
        changed.changes.push_back(GenotypeChange{.infection = infection, .locus = locus, .previous = previous});
        changed.stale = true;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::applyChanges(LikelihoodCache& current) {
        // Samplers that propose whole genomes set every locus, mostly to the value it already has, so only the loci
        // whose genotypes differ from those summed are changed loci.
        changedLoci_.clear();
        for (size_t c = 0; c < current.changes.size(); ++c) {
            const auto& change = current.changes[c];
            if (not child_->latentGenotype().contains(change.locus)) {
                continue;// terms are summed over the child's loci only
            }
            if (change.infection->latentGenotype(change.locus)->value() == change.previous) {
                continue;
            }
            if (std::ranges::find(changedLoci_, change.locus.get(), [](const auto& locus) { return locus.get(); }) == changedLoci_.end()) {
                changedLoci_.push_back(change.locus);
            }
        }
        // Swapping costs two locus terms per changed locus, recomputing one per locus of the child.
        const bool recompute = 2 * changedLoci_.size() > child_->loci().size();

        const size_t numClasses = current.classSize.size();
        std::array<const Genetics*, ParentSetMaxCardinality> previousCombination{};
        std::array<const Genetics*, ParentSetMaxCardinality> currentCombination{};
        bool childChanged = false;
        classChanged_.assign(numClasses, 0);

        // A class's genotypes are its member's. Only a member with recorded changes can differ from what was summed: a
        // class whose members all changed genome begins as the class they left (see regroup()).
        changedMembers_.clear();
        for (std::uint32_t k = 0; k < numClasses; ++k) {
            const auto* member = current.parents[current.classMember[k]];
            if (std::ranges::any_of(current.changes, [&](const auto& change) { return change.infection == member; })) {
                changedMembers_.emplace_back(member, k);
            }
        }

        // Calls update(position, combination, withLatent) for every entry involving a changed genotype.
        auto forEachAffected = [&](auto&& update) {
            if (childChanged) {
                update(0, std::span<const unsigned int>{}, true);
            }
            size_t position = 1;
            forEachCombination(numClasses, [&](const std::span<const unsigned int> combination) {
                const bool classesChanged = childChanged or std::any_of(combination.begin(), combination.end(), [&](const unsigned int k) { return classChanged_[k] != 0; });
                if (classesChanged) {
                    update(position, combination, false);
                    update(position + 1, combination, true);
                }
                position += 2;
            });
        };

        for (const auto& locus : changedLoci_) {
            // Every genotype at this locus as summed and as it is now.
            const Genetics currentChild = child_->latentGenotype(locus)->value();
            Genetics previousChild      = currentChild;
            if (not recompute) {
                childChanged = false;
                classChanged_.assign(numClasses, 0);
                previousClasses_.resize(numClasses);
                currentClasses_.resize(numClasses);
                for (size_t k = 0; k < numClasses; ++k) {
                    currentClasses_[k] = previousClasses_[k] = current.parents[current.classMember[k]]->latentGenotype(locus)->value();
                }
            }
            for (const auto& change : current.changes) {
                if (change.locus != locus) {
                    continue;
                }
                if (change.infection == child_.get()) {
                    previousChild = change.previous;
                    childChanged  = childChanged or previousChild != currentChild;
                } else {
                    for (const auto& [member, k] : changedMembers_) {
                        if (member == change.infection) {
                            if (not recompute) {
                                previousClasses_[k] = change.previous;
                            }
                            classChanged_[k] = classChanged_[k] or change.previous != member->latentGenotype(locus)->value();
                        }
                    }
                }
            }
            if (recompute) {
                continue;// entries are marked once every locus' changes are known
            }
            // The latent entries' tables are at the frequencies they were summed with.
            const FrequencyId frequencyId = current.frequencyIds.empty() ? 0 : current.frequencyIds[locusIndex_.at(locus.get())];

            forEachAffected([&](const size_t position, const std::span<const unsigned int> combination, const bool withLatent) {
                if (current.valid[position] == INVALID) {
                    return;// computed from scratch in value()
                }
                for (size_t t = 0; t < combination.size(); ++t) {
                    previousCombination[t] = &previousClasses_[combination[t]];
                    currentCombination[t]  = &currentClasses_[combination[t]];
                }
                const std::span<const Genetics* const> previousParents(previousCombination.data(), combination.size());
                const std::span<const Genetics* const> currentParents(currentCombination.data(), combination.size());
                if (withLatent) {
                    auto& sums = current.latentSums[position / 2];
                    if (sums.stoppedAt != nullptr) {
                        if (sums.stoppedAt == locus.get() and not memo_->locusTable(frequencyId, currentChild, currentParents).impossible) {
                            current.valid[position] = INVALID;
                        }
                        return;
                    }
                    sums.remove(memo_->locusTable(frequencyId, previousChild, previousParents));
                    sums.add(memo_->locusTable(frequencyId, currentChild, currentParents));
                } else {
                    auto& sums = current.sums[position];
                    if (sums.stoppedAt != nullptr) {
                        // Summed only up to a locus that makes the parent set impossible. It stays impossible unless that
                        // locus' term becomes possible, and is then computed from scratch.
                        if (sums.stoppedAt == locus.get() and not memo_->locusTerm(Variant::ObservedParents, currentChild, currentParents, nullptr).impossible) {
                            current.valid[position] = INVALID;
                        }
                        return;
                    }
                    sums.remove(memo_->locusTerm(Variant::ObservedParents, previousChild, previousParents, nullptr));
                    sums.add(memo_->locusTerm(Variant::ObservedParents, currentChild, currentParents, nullptr));
                }
                current.valid[position] = std::min(current.valid[position], SUMS_CHANGED);
                if (++current.updates[position] >= MAX_LOCUS_UPDATES) {
                    current.valid[position] = INVALID;
                }
            });
        }
        if (recompute) {
            const auto locusAddress = [](const LocusPtr& locus) { return locus.get(); };
            std::ranges::sort(changedLoci_, std::less{}, locusAddress);
            forEachAffected([&](const size_t position, std::span<const unsigned int>, const bool withLatent) {
                // An entry summed only up to a locus that makes it impossible stays impossible if that locus is unchanged.
                const auto* stoppedAt = withLatent ? static_cast<const void*>(current.latentSums[position / 2].stoppedAt) : static_cast<const void*>(current.sums[position].stoppedAt);
                if (current.valid[position] != INVALID and stoppedAt != nullptr and
                    not std::ranges::binary_search(changedLoci_, stoppedAt, std::less{}, [](const LocusPtr& locus) { return static_cast<const void*>(locus.get()); })) {
                    return;
                }
                current.valid[position] = INVALID;
            });
        }
        current.changes.clear();
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::syncFrequencyGeneration(LikelihoodCache& current) {
        if (current.frequencyGeneration == memo_->generation() and current.frequencyIds.size() == child_->loci().size()) {
            return;
        }
        current.frequencyIds.resize(child_->loci().size());
        for (size_t l = 0; l < child_->loci().size(); ++l) {
            current.frequencyIds[l] = memo_->frequencyId(child_->loci()[l]);
        }
        current.frequencyGeneration = memo_->generation();
        for (size_t position = 0; position < current.valid.size(); position += 2) {
            current.valid[position] = INVALID;
        }
        current.stale = true;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::applyFrequencyChanges(LikelihoodCache& current) {
        const size_t numClasses = current.classSize.size();
        std::array<const Genetics*, ParentSetMaxCardinality> combinationGenotypes{};
        for (const size_t l : pendingFrequencyLoci_) {
            const auto& locus           = child_->loci()[l];
            const FrequencyId updated   = memo_->frequencyId(locus);
            const FrequencyId previous  = current.frequencyIds[l];
            if (updated == previous) {
                continue;
            }
            current.frequencyIds[l] = updated;
            current.stale           = true;
            const Genetics& child   = child_->latentGenotype(locus)->value();
            currentClasses_.resize(numClasses);
            for (size_t k = 0; k < numClasses; ++k) {
                currentClasses_[k] = current.parents[current.classMember[k]]->latentGenotype(locus)->value();
            }
            auto swap = [&](const size_t position, const std::span<const unsigned int> combination) {
                if (current.valid[position] == INVALID) {
                    return;
                }
                for (size_t t = 0; t < combination.size(); ++t) {
                    combinationGenotypes[t] = &currentClasses_[combination[t]];
                }
                const std::span<const Genetics* const> parents(combinationGenotypes.data(), combination.size());
                auto& sums = current.latentSums[position / 2];
                if (sums.stoppedAt != nullptr) {
                    if (sums.stoppedAt == locus.get() and not memo_->locusTable(updated, child, parents).impossible) {
                        current.valid[position] = INVALID;
                    }
                    return;
                }
                sums.remove(memo_->locusTable(previous, child, parents));
                sums.add(memo_->locusTable(updated, child, parents));
                current.valid[position] = std::min(current.valid[position], SUMS_CHANGED);
                if (++current.updates[position] >= MAX_LOCUS_UPDATES) {
                    current.valid[position] = INVALID;
                }
            };
            swap(0, {});
            size_t position = 1;
            forEachCombination(numClasses, [&](const std::span<const unsigned int> combination) {
                swap(position + 1, combination);
                position += 2;
            });
        }
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    typename OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::LatentSums OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::latentSumsFromScratch(const LikelihoodCache& current, const std::span<const InfectionEventImpl* const> parents) {
        LatentSums sums{};
        const auto& loci = child_->loci();
        for (size_t l = 0; l < loci.size(); ++l) {
            const auto& locus = loci[l];
            for (size_t t = 0; t < parents.size(); ++t) {
                parentGenotypes_[t] = &parents[t]->latentGenotype(locus)->value();
            }
            const auto& table = memo_->locusTable(current.frequencyIds[l], child_->latentGenotype(locus)->value(), std::span<const Genetics* const>(parentGenotypes_.data(), parents.size()));
            sums.add(table);
            if (table.impossible) {
                sums.stoppedAt = locus.get();
                break;
            }
        }
        return sums;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::updateGroup(LikelihoodCache& current, LogSumExpGroup& group, const size_t first) {
        constexpr Likelihood negInf = -std::numeric_limits<Likelihood>::infinity();
        const size_t numEntries     = current.lliks.size();
        auto term                   = [&](const size_t i) { return current.valid[i] == INVALID ? negInf : current.lliks[i] + current.logMultiplicity[i]; };
        // When the entry holding the previous maximum is unchanged, every other unchanged entry is at most that
        // maximum, so only the newly computed entries can raise it.
        Likelihood maxTerm = negInf;
        size_t argmax      = NO_POSITION;
        if (group.weightArgmax < numEntries and current.valid[group.weightArgmax] == VALID) {
            assert(term(group.weightArgmax) == group.weightMax);
            maxTerm = group.weightMax;
            argmax  = group.weightArgmax;
            for (size_t i = first; i < numEntries; i += 2) {
                if (current.valid[i] == COMPUTED and term(i) > maxTerm) {
                    maxTerm = term(i);
                    argmax  = i;
                }
            }
        } else {
            for (size_t i = first; i < numEntries; i += 2) {
                if (const Likelihood t = term(i); argmax == NO_POSITION or t > maxTerm) {
                    maxTerm = t;
                    argmax  = i;
                }
            }
        }

        if (maxTerm == negInf) {
            group.logSum       = maxTerm;
            group.weightMax    = std::numeric_limits<Likelihood>::quiet_NaN();
            group.weightArgmax = NO_POSITION;
            return;
        }
        // Every weight is relative to the maximum, so a new maximum recomputes them all. Summed in position order.
        // Entries standing for no parent set weigh nothing.
        const bool rescale = maxTerm != group.weightMax;
        Likelihood sum     = 0;
        for (size_t i = first; i < numEntries; i += 2) {
            if (rescale or current.valid[i] != VALID) {
                current.weights[i] = std::exp(term(i) - maxTerm);
                if (current.valid[i] == COMPUTED) {
                    current.valid[i] = VALID;
                }
            }
            sum += current.weights[i];
        }
        group.weightMax    = maxTerm;
        group.weightArgmax = argmax;
        group.logSum       = maxTerm + std::log(sum);
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::combine(const LikelihoodCache& current) {
        // log( W0 + sum_S E_L[ W(S ∪ {L}) ] ): the latent parent is drawn from the Source whether or not it infected
        // the infection, so summing it out leaves the sets without it as they are.
        const Likelihood withoutLatent = current.groups[0].logSum;
        const Likelihood withLatent    = current.groups[1].logSum;
        const Likelihood maxTerm       = std::max(withoutLatent, withLatent);
        if (maxTerm == -std::numeric_limits<Likelihood>::infinity()) {
            return maxTerm;
        }
        return maxTerm + std::log(std::exp(withoutLatent - maxTerm) + std::exp(withLatent - maxTerm));
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    std::string OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::identifier() {
        auto out = fmt::format("OrderBasedTransmissionProcessV3<{}>", child_->id());
        return out;
    }


    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::value() {
        if (this->isDirty() and not cache().stale and not parentsMayHaveChanged_ and not frequenciesMayHaveChanged_) {
            // Nothing the entries depend on changed: recombine the groups' sums.
            this->value_ = combine(cache());
            this->setClean();
        }
        if (this->isDirty()) {
            // Genome ids, used to group parents and compute entries, stay valid until the next batch.
            memo_->beginBatch();
            // Read every evaluation, so that a change to it is never accepted unread (see the listener).
            const auto coiLogProbs = coiProb_->value();
            const auto& ps           = parentSet_->members();
            const bool parentChanged = std::ranges::any_of(cache().changes, [&](const auto& change) {
                return change.infection != child_.get();
            });
            if (parentsMayHaveChanged_ or parentChanged or cache().valid.empty()) {
                regroup(ps, parentsMayHaveChanged_);
            }
            parentsMayHaveChanged_ = false;
            auto& current          = mutableCache();

            if (null_model_) {
                current.changes.clear();
                for (size_t position = 0; position < current.valid.size(); ++position) {
                    if (current.logMultiplicity[position] != -std::numeric_limits<Likelihood>::infinity()) {
                        current.lliks[position] = 0.0;
                        current.valid[position] = COMPUTED;
                    }
                }
            } else {
                syncFrequencyGeneration(current);
                applyChanges(current);
                if (frequenciesMayHaveChanged_) {
                    applyFrequencyChanges(current);
                }

                // Genome ids are resolved once per evaluation rather than per parent set; parents' lazily, since
                // usually only a few entries need computing from scratch.
                typename MemoImpl::ResolvedGenomes resolved{};
                bool childResolved = false;
                auto resolveChild  = [&]() {
                    if (not childResolved) {
                        resolved.childLocusOrder = memo_->locusOrderOf(child_);
                        resolved.child           = memo_->genomeOf(child_);
                        childResolved            = true;
                    }
                };
                bool genomesReset = false;
                auto parentGenome = [&](const unsigned int idx) {
                    if (not genomesReset) {
                        parentGenomes_.assign(ps.size(), std::nullopt);
                        genomesReset = true;
                    }
                    if (not parentGenomes_[idx]) {
                        parentGenomes_[idx] = memo_->genomeOf(ps.begin()[idx]);
                    }
                    return *parentGenomes_[idx];
                };
                auto setSums = [&](const size_t position, const StrainSums& sums) {
                    current.sums[position]    = sums;
                    current.updates[position] = 0;
                    current.valid[position]   = SUMS_CHANGED;
                };
                auto setLatentSums = [&](const size_t position, const LatentSums& sums) {
                    current.latentSums[position / 2] = sums;
                    current.updates[position]        = 0;
                    current.valid[position]          = SUMS_CHANGED;
                };

                // The outer term: the parent-set-size prior. psp holds probabilities, not log-probabilities.
                const auto parentSetSizeProbs = psp_->value();
                std::array<Likelihood, ParentSetMaxCardinality + 2> logParentSetSizeProb{};
                for (size_t k = 1; k < logParentSetSizeProb.size(); ++k) {
                    logParentSetSizeProb[k] = std::log(parentSetSizeProbs(k));
                }
                // Brings an entry from any state with current sums up to COMPUTED, in the same operation order as
                // computing it from scratch.
                auto finish = [&](const size_t position, const size_t numParents) {
                    if (current.valid[position] == SUMS_CHANGED) {
                        current.summed[position] = position % 2 == 0 ? memo_->latentSummedLogLikelihood(current.latentSums[position / 2], numParents, coiLogProbs)
                                                                     : memo_->strainSummedLogLikelihood(current.sums[position], numParents);
                        current.valid[position]  = OUTER_CHANGED;
                    }
                    if (current.valid[position] == OUTER_CHANGED) {
                        Likelihood llik = current.summed[position];
                        if (llik != -std::numeric_limits<Likelihood>::infinity()) {
                            llik += logParentSetSizeProb[numParents];
                        }
                        current.lliks[position] = llik;
                        current.valid[position] = COMPUTED;
                    }
                };

                // Calculate the single latent parent case
                if (current.valid[0] == INVALID) {
                    resolveChild();
                    setLatentSums(0, latentSumsFromScratch(current, {}));
                }
                finish(0, 1);

                // Every combination of classes, computing from scratch only the sums not already cached. Entries
                // standing for no parent set (a pair within a class of one) are left out.
                auto needed = [&](const size_t position) {
                    return current.valid[position] == INVALID and current.logMultiplicity[position] != -std::numeric_limits<Likelihood>::infinity();
                };
                ParentSet tmpPs_{};
                size_t position = 1;
                forEachCombination(current.classSize.size(), [&](const std::span<const unsigned int> combination) {
                    if (needed(position) or needed(position + 1)) {
                        // Members standing in for the classes, as indices into the parent set: two of one class for a
                        // pair within it.
                        std::array<unsigned int, 2> members{current.classMember[combination[0]], 0};
                        if (combination.size() == 2) {
                            if (combination[0] != combination[1]) {
                                members[1] = current.classMember[combination[1]];
                            } else {
                                for (unsigned int i = 0; i < current.parents.size(); ++i) {
                                    if (current.parentClass[i] == combination[0] and i != members[0]) {
                                        members[1] = i;
                                        break;
                                    }
                                }
                            }
                            // Ascending indices into the uid-ordered parent set, so this matches the order of resolved.parents.
                            if (members[0] > members[1]) {
                                std::swap(members[0], members[1]);
                            }
                        }
                        resolveChild();
                        tmpPs_.clear();
                        for (size_t t = 0; t < combination.size(); ++t) {
                            resolved.parents[t] = parentGenome(members[t]);
                            tmpPs_.insert(ps.begin()[members[t]]);
                        }
                        if (needed(position)) {
                            setSums(position, memo_->parentSetSums(Variant::ObservedParents, child_, &tmpPs_, nullptr, resolved));
                        }
                        if (needed(position + 1)) {
                            std::array<const InfectionEventImpl*, 2> chosen{};
                            for (size_t t = 0; t < combination.size(); ++t) {
                                chosen[t] = ps.begin()[members[t]].get();
                            }
                            setLatentSums(position + 1, latentSumsFromScratch(current, std::span<const InfectionEventImpl* const>(chosen.data(), combination.size())));
                        }
                    }
                    finish(position, combination.size());
                    finish(position + 1, combination.size() + 1);
                    position += 2;
                });
            }

            updateGroup(current, current.groups[0], 1);
            updateGroup(current, current.groups[1], 0);
            current.stale = false;
            this->value_  = combine(current);
            frequenciesMayHaveChanged_ = false;
            for (const size_t l : pendingFrequencyLoci_) {
                frequencyPending_[l] = 0;
            }
            pendingFrequencyLoci_.clear();
            assert(this->value_ < std::numeric_limits<Likelihood>::infinity());

            this->setClean();
        }
        if (null_model_) {
            return 0;
        }
        return this->value_;
    }


    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::postSaveState([[maybe_unused]] int savedStateId) {
        slots_[cacheIndex_ + 1] = slots_[cacheIndex_];
        cacheIndex_++;
    }


    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::postRestoreState([[maybe_unused]] int savedStateId) {
        cacheIndex_--;
    }


    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::postAcceptState() {
        if (slots_[cacheIndex_] != 0) {
            std::swap(caches_[0], caches_[slots_[cacheIndex_]]);
        }
        cacheIndex_ = 0;
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::peek() noexcept {
        if (null_model_) {
            return 0;
        }
        return Computation::peek();
    }

    template<int ParentSetMaxCardinality, typename MemoImpl, typename COIProbabilityImpl, typename ParentSetSizeLikelihoodImpl, typename InfectionEventImpl, typename ParentSetImpl>
    ParentSetDist<InfectionEventImpl> OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, MemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEventImpl, ParentSetImpl>::calcParentSetDist() {
        ParentSetDist<InfectionEventImpl> dist{};
        dist.totalLlik          = this->value();
        const auto& current     = cache();
        const auto tmpPs        = parentSet_->value();
        const size_t totalNodes = tmpPs.size();
        ParentSet ps{};
        ps.insert(latentParent_);
        dist.parentSetLliks.push_back(std::make_pair(current.lliks[0], ps));

        // Each parent set has the likelihood of its classes' entry.
        std::array<unsigned int, 2> classes{};
        auto addParentSet = [&](const std::span<const unsigned int> parents) {
            ps.clear();
            for (size_t t = 0; t < parents.size(); ++t) {
                ps.insert(tmpPs.begin()[parents[t]]);
                classes[t] = current.parentClass[parents[t]];
            }
            if (parents.size() == 2 and classes[0] > classes[1]) {
                std::swap(classes[0], classes[1]);
            }
            const size_t position = 1 + 2 * combinationIndex(std::span<const unsigned int>(classes.data(), parents.size()));
            dist.parentSetLliks.push_back(std::make_pair(current.lliks[position], ps));
            ps.insert(latentParent_);
            dist.parentSetLliks.push_back(std::make_pair(current.lliks[position + 1], ps));
        };
        // In the order of the parent sets themselves: all single parents, then all pairs.
        std::array<unsigned int, 2> parents{};
        for (unsigned int i = 0; i < totalNodes; ++i) {
            parents[0] = i;
            addParentSet(std::span<const unsigned int>(parents.data(), 1));
        }
        if constexpr (ParentSetMaxCardinality == 2) {
            for (unsigned int i = 0; i + 1 < totalNodes; ++i) {
                parents[0] = i;
                for (unsigned int j = i + 1; j < totalNodes; ++j) {
                    parents[1] = j;
                    addParentSet(std::span<const unsigned int>(parents.data(), 2));
                }
            }
        }

        if (null_model_) {
            Likelihood totalLlik = 0;
            for (const auto& [llik, ps] : dist.parentSetLliks) {
                totalLlik += std::exp(llik);
            }
            dist.totalLlik = std::log(totalLlik);
        }

        return dist;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H
