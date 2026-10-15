//
// The Ordering: every Infection in the order it was acquired, derived from observation time minus infection duration,
// and each Infection's Candidate parents: the parents the study permits it that precede it in the order. A duration
// change moves its Infection one neighbour at a time until the order is sorted again; each crossing tells the one or
// two Infections whose Candidate parents gained or lost a member (the keyed candidates_changed event), and the view
// candidateParents() returns is recomputed lazily, once per change, for the Infections that read it. A genotype change
// of an Infection is relayed to each Infection it is a Candidate parent of (the keyed candidate_genotype_changed event).
// The order is saved and restored with the durations; a restore or a checkpoint's restoreOrder rebuilds the positions.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Dependency.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/containers/ParentSet.h"
#include "core/containers/UidLess.h"

#include <boost/container/flat_set.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace transmission_nets::core::computation {

    template<typename InfectionEventImpl>
    class ObservationTimeDerivedOrdering : public Computation<std::vector<std::shared_ptr<InfectionEventImpl>>>,
                                           public abstract::Observable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>,
                                           public abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<std::shared_ptr<InfectionEventImpl>>> {
        using Infection = std::shared_ptr<InfectionEventImpl>;
        using ChangedCallback = std::function<void()>;
        // Fired for an Infection whose Candidate parents gained or lost a member, as the order moves.
        CREATE_KEYED_EVENT(candidates_changed, Infection, ChangedCallback);
        // Fired for an Infection one of whose Candidate parents changed genotype.
        CREATE_KEYED_EVENT(candidate_genotype_changed, Infection, ChangedCallback);

    public:
        using ParentSet        = containers::ParentSet<InfectionEventImpl>;
        using PermittedParents = std::function<std::vector<Infection>(const Infection&)>;

        ObservationTimeDerivedOrdering() noexcept;

        // Adds `infections`, in infection-time order, with `permittedParents(infection)` the parents the study
        // permits each; only those can be Candidate parents.
        void addElements(const std::vector<Infection>& infections, const PermittedParents& permittedParents);

        std::vector<Infection> value() override;

        // Puts the elements in `order`, which holds exactly the current elements, as when restoring a checkpoint:
        // Infections with equal infection times keep their saved order. Notifies no listener, so it is for use
        // before anything depends on the ordering.
        void restoreOrder(std::vector<Infection> order) noexcept {
            assert(order.size() == this->value_.size());
            this->value_ = std::move(order);
            orderChanged();
            this->setDirty();
        }

        // `child`'s Candidate parents, in uid order: the parents the study permits it that precede it. The reference
        // is valid until the order next changes.
        const ParentSet& candidateParents(const Infection& child) {
            auto& cached = candidates_[uid(child)];
            if (cached.version != version_) {
                cached.parents.clear();
                for (const auto& parent : permitted_[uid(child)]) {
                    if (positions_[uid(parent)] < positions_[uid(child)]) {
                        cached.parents.push_back(parent);
                    }
                }
                cached.set     = ParentSet(boost::container::ordered_unique_range, cached.parents.begin(), cached.parents.end());
                cached.version = version_;
            }
            return cached.set;
        }

        // The parents the study permits `child`, in uid order, whether or not they precede it.
        const std::vector<Infection>& permittedParents(const Infection& child) const {
            return permitted_[uid(child)];
        }

        // Advances with every change to the order, restores included; a Candidate-parent view read at one version is
        // current until the next.
        [[nodiscard]] std::uint64_t version() const noexcept { return version_; }

        // Where `infection` stands, counting from the first acquired.
        [[nodiscard]] std::size_t position(const Infection& infection) const { return positions_[uid(infection)]; }

    protected:
        friend class abstract::Cacheable<ObservationTimeDerivedOrdering<InfectionEventImpl>>;
        friend class abstract::Checkpointable<ObservationTimeDerivedOrdering<InfectionEventImpl>, std::vector<Infection>>;

    private:
        struct Candidates {
            std::uint64_t version = std::numeric_limits<std::uint64_t>::max();
            std::vector<Infection> parents;// scratch, in uid order
            ParentSet set;
        };

        void infectionDurationChanged(const Infection& ref);
        void genotypeChanged(const Infection& parent);
        void orderChanged();
        void crossed(const Infection& nowBefore, const Infection& nowAfter);
        bool permits(const Infection& child, const Infection& parent) const;

        static std::size_t uid(const Infection& infection) { return static_cast<std::size_t>(infection->uid()); }

        template<typename Slot>
        Slot& byUid(std::vector<Slot>& slots, const Infection& infection) {
            if (slots.size() <= uid(infection)) {
                slots.resize(uid(infection) + 1);
            }
            return slots[uid(infection)];
        }

        std::vector<std::size_t> positions_{};           // by uid
        std::vector<std::vector<Infection>> permitted_{};// by child's uid, in uid order
        std::vector<std::vector<Infection>> children_{}; // by parent's uid: the Infections it is permitted to, in input order
        std::vector<Candidates> candidates_{};           // by child's uid
        std::uint64_t version_ = 0;
    };

    template<typename InfectionEventImpl>
    ObservationTimeDerivedOrdering<InfectionEventImpl>::ObservationTimeDerivedOrdering() noexcept {
        this->addPostRestoreHook([this](int) { orderChanged(); });
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::addElements(const std::vector<Infection>& infections, const PermittedParents& permittedParents) {
        for (const auto& ref : infections) {
            this->value_.push_back(ref);
            register_candidates_changed_listener_key(ref);
            register_candidate_genotype_changed_listener_key(ref);
            byUid(positions_, ref);
            byUid(candidates_, ref);
            byUid(children_, ref);
            auto& permitted = byUid(permitted_, ref);
            permitted       = permittedParents(ref);
            std::ranges::sort(permitted, containers::UidLess{});
            abstract::dependsOn(this, ref->infectionDuration(), [this, ref]() { this->infectionDurationChanged(ref); });
        }
        for (const auto& child : infections) {
            for (const auto& parent : permitted_[uid(child)]) {
                byUid(children_, parent).push_back(child);
            }
        }
        for (const auto& ref : infections) {
            // A latent genotype change leaves the ordering as it is, so the ordering is not saved with it; it is
            // relayed to the Infections `ref` is a Candidate parent of.
            ref->add_post_change_listener([this, ref]() { this->genotypeChanged(ref); });
        }

        std::sort(this->value_.begin(),
                  this->value_.end(),
                  [](const Infection& a, const Infection& b) { return a->infectionTime() < b->infectionTime(); });
        orderChanged();
        this->setDirty();
    }

    template<typename InfectionEventImpl>
    std::vector<std::shared_ptr<InfectionEventImpl>> ObservationTimeDerivedOrdering<InfectionEventImpl>::value() {
        this->setClean();
        return this->value_;
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::genotypeChanged(const Infection& parent) {
        for (const auto& child : children_[uid(parent)]) {
            if (positions_[uid(parent)] < positions_[uid(child)]) {
                keyed_notify_candidate_genotype_changed(child);
            }
        }
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::orderChanged() {
        for (std::size_t at = 0; at < this->value_.size(); ++at) {
            positions_[uid(this->value_[at])] = at;
        }
        ++version_;
    }

    template<typename InfectionEventImpl>
    bool ObservationTimeDerivedOrdering<InfectionEventImpl>::permits(const Infection& child, const Infection& parent) const {
        const auto& permitted = permitted_[uid(child)];
        return std::binary_search(permitted.begin(), permitted.end(), parent, containers::UidLess{});
    }

    // `nowBefore` has just crossed `nowAfter`, so it precedes it from here on: if the study permits it, `nowAfter`
    // gains it as a Candidate parent, and `nowBefore` loses `nowAfter`.
    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::crossed(const Infection& nowBefore, const Infection& nowAfter) {
        if (permits(nowAfter, nowBefore)) {
            keyed_notify_candidates_changed(nowAfter);
        }
        if (permits(nowBefore, nowAfter)) {
            keyed_notify_candidates_changed(nowBefore);
        }
    }

    template<typename InfectionEventImpl>
    void ObservationTimeDerivedOrdering<InfectionEventImpl>::infectionDurationChanged(const Infection& ref) {
        double refInfectionTime = ref->infectionTime();

        // Find the index of the reference in the vector.
        size_t refIdx = positions_[uid(ref)];

        // Find the infection time of the infection just after the reference.
        // If the reference is the last element, then the infection time is a high number.
        double rightTime = std::numeric_limits<double>::max();
        if (refIdx < this->value_.size() - 1) {
            rightTime = this->value_[refIdx + 1]->infectionTime();
        }

        // Find the infection time of the infection just before the reference.
        // If the reference is the first element, then the infection time is a low number.
        double leftTime = std::numeric_limits<double>::min();
        if (refIdx > 0) {
            leftTime = (this->value_[refIdx - 1])->infectionTime();
        }

        // If the reference infection time is less than the left time or greater than the right time, then the reference is
        // out of order.
        bool unsorted = refInfectionTime < leftTime or refInfectionTime > rightTime;

        if (unsorted) {
            this->setDirty();
            ++version_;
            // if the reference infection time is less than the left time, then the reference is moved to the left.
            if (refInfectionTime < leftTime) {
                while (unsorted) {
                    // Make sure the reference is not the first element.
                    if (refIdx > 0) {
                        // Move the reference to the left.
                        auto& curr = this->value_[refIdx];
                        auto& prev = this->value_[refIdx - 1];
                        std::swap(curr, prev);
                        positions_[uid(this->value_[refIdx])]     = refIdx;
                        positions_[uid(this->value_[refIdx - 1])] = refIdx - 1;

                        // The reference now precedes the element it passed.
                        crossed(this->value_[refIdx - 1], this->value_[refIdx]);

                        // Update the reference index.
                        --refIdx;

                        // check if the reference is sorted.
                        if (refIdx == 0 or this->value_[refIdx - 1]->infectionTime() < refInfectionTime) {
                            unsorted = false;
                        }
                    } else {
                        unsorted = false;
                    }
                }
            } else {
                while (unsorted) {
                    // Make sure the reference is not the last element.
                    if (refIdx < this->value_.size() - 1) {
                        // Move the reference to the right.
                        auto& curr = this->value_[refIdx];
                        auto& next = this->value_[refIdx + 1];
                        std::swap(curr, next);
                        positions_[uid(this->value_[refIdx])]     = refIdx;
                        positions_[uid(this->value_[refIdx + 1])] = refIdx + 1;

                        // The element the reference passed now precedes it.
                        crossed(this->value_[refIdx], this->value_[refIdx + 1]);

                        // Update the reference index.
                        ++refIdx;

                        // check if the reference is sorted.
                        if (refIdx == this->value_.size() - 1 or this->value_[refIdx + 1]->infectionTime() > refInfectionTime) {
                            unsorted = false;
                        }
                    } else {
                        unsorted = false;
                    }
                }
            }
        }
    }

}// namespace transmission_nets::core::computation

#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONTIMEDERIVEDORDERING_H
