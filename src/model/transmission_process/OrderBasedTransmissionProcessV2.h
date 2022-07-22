//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV2_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV2_H

#include "core/computation/OrderDerivedParentSet.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/io/serialize.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"


#include <boost/container/flat_set.hpp>
#include <boost/range/adaptors.hpp>
#include <fmt/core.h>
#include <fmt/ranges.h>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<typename InfectionEventImpl>
    struct ParentSetDist {
        Likelihood sourceLlik = 0;
        std::vector<std::pair<Likelihood, core::containers::ParentSet<InfectionEventImpl>>> parentSetLliks{};
        Likelihood totalLlik = 0;
    };

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    class OrderBasedTransmissionProcessV2 : public core::computation::PartialLikelihood {

        /*
         * The order based transmission process considers the set of all possible parent sets under the given ordering.
         * Assumes each infection has a potential latent unobserved parent that is drawn from the background distribution.
         */

        using ListenerIdMap = boost::container::flat_map<std::shared_ptr<InfectionEventImpl>, core::abstract::ListenerId_t>;

    public:
        //// Node Transmission Process changes -> recalculate completely
        //// Child updated -> recalculate completely
        //// Parent Set adds node -> add node combo to update set, register listener
        //// Parent Set removes node -> subtract node combos immediately, remove listener
        //// Node updated in parent set -> subtract node combos immediately, add node combo to update set
        OrderBasedTransmissionProcessV2(std::shared_ptr<NodeTransmissionProcessImpl> ntp, std::shared_ptr<InfectionEventImpl> child, std::shared_ptr<ParentSetImpl> parent_set, std::shared_ptr<InfectionEventImpl> latent_parent);

        Likelihood value() override;

        Likelihood recalculate(bool verbose = false);

        Likelihood calculateParentLogLikelihoodContribution(std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others);

        Likelihood calculateParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps);

        ParentSetDist<InfectionEventImpl> calcParentSetDist();

        std::string identifier() override;

        std::shared_ptr<NodeTransmissionProcessImpl> ntp_;
        std::shared_ptr<InfectionEventImpl> child_;
        std::shared_ptr<ParentSetImpl> parentSet_;
        std::shared_ptr<InfectionEventImpl> latentParent_;

    private:
        void nodeTransmissionProcessSetDirty();

        void childSetDirty();

        void addParent(std::shared_ptr<InfectionEventImpl> parent);

        void addParentListeners(std::shared_ptr<InfectionEventImpl> parent);

        void removeParent(std::shared_ptr<InfectionEventImpl> parent);

        void removeParentListeners(std::shared_ptr<InfectionEventImpl> parent);

        void preParentUpdated(std::shared_ptr<InfectionEventImpl> parent);

        void postParentUpdated(std::shared_ptr<InfectionEventImpl> parent);

        void postSaveState();

        void postAcceptState();

        void postRestoreState();

        ListenerIdMap preChangeListenerIdMap{};
        ListenerIdMap postChangeListenerIdMap{};
        ListenerIdMap saveStateListenerIdMap{};
        ListenerIdMap acceptStateListenerIdMap{};
        ListenerIdMap restoreStateListenerIdMap{};

        using InfectionEventSet = boost::container::flat_set<std::shared_ptr<InfectionEventImpl>>;

        // Track parent deltas between save and accept/restore
        InfectionEventSet addedParents_{};
        InfectionEventSet removedParents_{};

        std::vector<InfectionEventSet> addedParentsCache_{};
        std::vector<InfectionEventSet> removedParentsCache_{};

        InfectionEventSet calculated_{};
        InfectionEventSet toCalculate_{};
        std::vector<Likelihood> toSubtract_;

        std::deque<InfectionEventSet> calculatedCache_{};
        std::deque<InfectionEventSet> toCalculateCache_{};
        std::deque<std::vector<Likelihood>> toSubtractCache_{};


        // helper vars for calculating parent likelihood contributions
        core::containers::ParentSet<InfectionEventImpl> tmpPs_{};
        std::vector<Likelihood> parentLikelihoodContribution_{};
        core::utils::generators::CombinationIndicesGenerator cs_;
        std::vector<Likelihood> toAdd_{};

    public:
        Likelihood peek() noexcept override;
    };


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::OrderBasedTransmissionProcessV2(
            std::shared_ptr<NodeTransmissionProcessImpl> ntp, std::shared_ptr<InfectionEventImpl> child, std::shared_ptr<ParentSetImpl> parent_set, std::shared_ptr<InfectionEventImpl> latent_parent) : ntp_(std::move(ntp)), child_(std::move(child)), parentSet_(std::move(parent_set)), latentParent_(std::move(latent_parent)) {

        ntp_->add_set_dirty_listener([=, this]() { nodeTransmissionProcessSetDirty(); });
        ntp_->registerCacheableCheckpointTarget(this);

        child_->add_post_change_listener([=, this]() { childSetDirty(); });
        child_->registerCacheableCheckpointTarget(this);

        parentSet_->add_element_added_listener([=, this](std::shared_ptr<InfectionEventImpl> parent) { addParent(parent); });
        parentSet_->add_element_removed_listener([=, this](std::shared_ptr<InfectionEventImpl> parent) { removeParent(parent); });
        parentSet_->registerCacheableCheckpointTarget(this);

        this->addPostSaveHook([=, this]() { this->postSaveState(); });
        this->addPostRestoreHook([=, this]() { this->postRestoreState(); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        // Register the latent parent and the parents in the derived parent set
        addParent(latentParent_);
        for (auto& parent : parentSet_->value()) {
            addParent(parent);
        }

        // Clear the added parents tracker on first initialization
        addedParents_.clear();
        this->value_ = -std::numeric_limits<Likelihood>::infinity();
        this->setDirty();
        this->value();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    std::string OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::identifier() {
        return "OrderBasedTransmissionProcessV2";
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::value() {
        if (this->isDirty()) {
            toAdd_.clear();
            Likelihood maxLlik = -std::numeric_limits<Likelihood>::infinity();
            Likelihood totalToSubtract = 0.0;
            if (!toSubtract_.empty()) {
                totalToSubtract = core::utils::logSumExp(toSubtract_);
            }

            if (std::abs(this->value_ - totalToSubtract) < 1e-24) {
                // numerical precision issues when subtracting two very close numbers
                this->value_ = recalculate();
            } else {
                if (totalToSubtract != 0.0) {
                    this->value_ = core::utils::logDiffExp(this->value_, totalToSubtract);
                }

                toAdd_.push_back(this->value_);
                maxLlik = std::max(maxLlik, toAdd_.back());
                for (const auto& parent : toCalculate_) {
                    assert(!(calculated_.contains(parent)));
                    toAdd_.push_back(calculateParentLogLikelihoodContribution(parent, calculated_));
                    maxLlik = std::max(maxLlik, toAdd_.back());
                    calculated_.insert(parent);
                }


                this->value_ = core::utils::logSumExpKnownMax(toAdd_.begin(), toAdd_.end(), maxLlik);

                assert(this->value_ < std::numeric_limits<Likelihood>::infinity());

#ifdef DEBUG_LIKELIHOOD

                auto tmp = recalculate();
                if (this->value_ > -std::numeric_limits<Likelihood>::infinity() and (std::abs(tmp - this->value_) > 1)) {
                    fmt::print("Err: {}, {}, {}\n", std::abs(tmp - this->value_), tmp, this->value_);
                    fmt::print("{} Likelihood mismatch OBTP: {}, {}\n", fmt::ptr(this), tmp, this->value_);
                    auto second_validate = recalculate(true);
                    fmt::print("Second Validation: {}\n", second_validate);
                    fmt::print("Child Node: {}\n", child_->id());
                    fmt::print("Total Parents: {}\n", parentSet_->value().size());
                    fmt::print("stp Value: {}\n", stp_->value());
                    fmt::print("To Add: {}\n", fmt::join(toAdd_, ", "));
                    fmt::print("To Subtract: {}\n", fmt::join(toSubtract_, ", "));
                    fmt::print("Total to subtract: {0:.36f}\n", totalToSubtract);
                    std::vector<std::string> removedParentIDs{};
                    std::vector<std::string> currentParentIDs{};
                    std::vector<std::string> toCalculateIDs{};
                    std::vector<std::string> calculatedParents{};
                    std::transform(removedParents_.begin(), removedParents_.end(), std::back_inserter(removedParentIDs), [](auto p) { return p->id(); });
                    std::transform(toCalculate_.begin(), toCalculate_.end(), std::back_inserter(toCalculateIDs), [](auto p) { return p->id(); });
                    std::transform(
                            calculated_.begin(), calculated_.end(), std::back_inserter(calculatedParents), [](auto p) { return p->id(); });
                    for (const auto& parent : parentSet_->value()) {
                        currentParentIDs.push_back(parent->id());
                    }
                    fmt::print("Removed Parents: {}\n", fmt::join(removedParentIDs, ", "));
                    fmt::print("Current Parents: {}\n", fmt::join(currentParentIDs, ", "));
                    fmt::print("To Calculate: {}\n", fmt::join(toCalculateIDs, ", "));
                    fmt::print("Calculated Parents: {}\n", fmt::join(calculatedParents, ", "));

                    exit(1);
                }
#endif
            }
            this->setClean();
            toSubtract_.clear();
            toCalculate_.clear();
        }

        return this->value_;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::recalculate(bool verbose) {
        /*
         * function that completely recalculates the likelihood.
         */
        auto val = -std::numeric_limits<Likelihood>::infinity();
        InfectionEventSet tmpCalculated_{};
        std::vector<Likelihood> tmpToAdd_{};
        Likelihood maxLlik = val;

        std::vector<std::string> parentIDs{};

        tmpToAdd_.push_back(calculateParentLogLikelihoodContribution(latentParent_, tmpCalculated_));
        tmpCalculated_.insert(latentParent_);
        maxLlik = std::max(maxLlik, tmpToAdd_.back());

        for (const auto& parent : parentSet_->value()) {
            tmpToAdd_.push_back(calculateParentLogLikelihoodContribution(parent, tmpCalculated_));
            maxLlik = std::max(maxLlik, tmpToAdd_.back());
            tmpCalculated_.insert(parent);
            parentIDs.push_back(parent->id());
        }
        maxLlik = std::max(maxLlik, tmpToAdd_.back());

        if (verbose) {
            fmt::print("Validate to Add: {}\n", fmt::join(tmpToAdd_, ", "));
            fmt::print("Validation Current Parents: {}\n", fmt::join(parentIDs, ", "));
        }

        val         = core::utils::logSumExpKnownMax(tmpToAdd_.begin(), tmpToAdd_.end(), maxLlik);
        calculated_ = std::move(tmpCalculated_);
        return val;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::nodeTransmissionProcessSetDirty() {
        /*
         * Completely recalculate the likelihood
         */

        if (!this->isDirty()) {
            this->setDirty();
        }

        toCalculate_.clear();
        toSubtract_.clear();
        calculated_.clear();

        toCalculate_.insert(latentParent_);
        for (auto& p : parentSet_->value()) {
            toCalculate_.insert(p);
        }

        this->value_ = -std::numeric_limits<Likelihood>::infinity();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::childSetDirty() {
        /*
         * Completely recalculate the likelihood
         */

        if (!this->isDirty()) {
            this->setDirty();
        }

        toCalculate_.clear();
        toSubtract_.clear();
        calculated_.clear();

        toCalculate_.insert(latentParent_);
        for (auto& p : parentSet_->value()) {
            toCalculate_.insert(p);
        };
        this->value_ = -std::numeric_limits<Likelihood>::infinity();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::addParent(
            std::shared_ptr<InfectionEventImpl> parent) {
        /*
         * Mark new parent to be calculated, add listeners to parent
         */

        addParentListeners(parent);
        addedParents_.insert(parent);
        toCalculate_.insert(parent);

        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::addParentListeners(
            std::shared_ptr<InfectionEventImpl> parent) {
        const auto preChangeListenerId    = parent->add_pre_change_listener([=, this]() { preParentUpdated(parent); });
        const auto postChangeListenerId   = parent->add_post_change_listener([=, this]() { postParentUpdated(parent); });
        const auto saveStateListenerId    = parent->add_save_state_listener([=, this](const std::string& savedStateId) { saveState(savedStateId); });
        const auto acceptStateListenerId  = parent->add_accept_state_listener([=, this]() { acceptState(); });
        const auto restoreStateListenerId = parent->add_restore_state_listener([=, this](const std::string& savedStateId) { restoreState(savedStateId); });

        assert(!preChangeListenerIdMap.contains(parent));
        preChangeListenerIdMap[parent]    = preChangeListenerId;
        postChangeListenerIdMap[parent]   = postChangeListenerId;
        saveStateListenerIdMap[parent]    = saveStateListenerId;
        acceptStateListenerIdMap[parent]  = acceptStateListenerId;
        restoreStateListenerIdMap[parent] = restoreStateListenerId;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::removeParent(
            std::shared_ptr<InfectionEventImpl> parent) {

        calculated_.erase(parent);

        toSubtract_.push_back(calculateParentLogLikelihoodContribution(parent, calculated_));

        removeParentListeners(parent);
        removedParents_.insert(parent);

        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::removeParentListeners(
            std::shared_ptr<InfectionEventImpl> parent) {
        parent->remove_pre_change_listener(preChangeListenerIdMap.at(parent));
        parent->remove_post_change_listener(postChangeListenerIdMap.at(parent));
        parent->remove_save_state_listener(saveStateListenerIdMap.at(parent));
        parent->remove_accept_state_listener(acceptStateListenerIdMap.at(parent));
        parent->remove_restore_state_listener(restoreStateListenerIdMap.at(parent));
        preChangeListenerIdMap.erase(parent);
        postChangeListenerIdMap.erase(parent);
        saveStateListenerIdMap.erase(parent);
        acceptStateListenerIdMap.erase(parent);
        restoreStateListenerIdMap.erase(parent);
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::preParentUpdated(
            std::shared_ptr<InfectionEventImpl> parent) {
        // Before the parent is updated, we need to remove the contribution
        // of the parent from the likelihood.
        if (calculated_.contains(parent)) {
            calculated_.erase(parent);
            toSubtract_.push_back(calculateParentLogLikelihoodContribution(parent, calculated_));
            toCalculate_.insert(parent);
        }
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::postParentUpdated(
            [[maybe_unused]] std::shared_ptr<InfectionEventImpl> parent) {
        // After the parent is updated, set the likelihood to dirty.
        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::calculateParentLogLikelihoodContribution(
            std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others) {
        const int otherNodesSize = others.size();
        parentLikelihoodContribution_.clear();

        // Calculate the single parent case
        tmpPs_.clear();
        tmpPs_.insert(parent);
        parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, tmpPs_));
        Likelihood max_llik = parentLikelihoodContribution_.back();
        for (int i = 1; i < ParentSetMaxCardinality and i <= otherNodesSize; ++i) {
            cs_.reset(otherNodesSize, i);
            while (!cs_.completed) {
                tmpPs_.clear();
                for (const auto& idx : cs_.curr) {
                    tmpPs_.insert(others.begin()[idx]);
                }
                tmpPs_.insert(parent);
                parentLikelihoodContribution_.push_back(calculateParentSetLogLikelihood(tmpPs_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                cs_.next();
            }
        }

        Likelihood val = core::utils::logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
        assert(val < std::numeric_limits<Likelihood>::infinity());
        return val;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::postSaveState() {
        calculatedCache_.emplace_back(calculated_);
        toCalculateCache_.emplace_back(toCalculate_);
        toSubtractCache_.emplace_back(toSubtract_);

        addedParentsCache_.emplace_back(addedParents_);
        removedParentsCache_.emplace_back(removedParents_);

        addedParents_.clear();
        removedParents_.clear();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::postRestoreState() {
        for (const auto& parent : addedParents_) {
            removeParentListeners(parent);
        }

        for (const auto& parent : removedParents_) {
            addParentListeners(parent);
        }

        assert(!(calculatedCache_.empty()));
        assert(!(removedParentsCache_.empty()));
        assert(!(addedParentsCache_.empty()));

        calculated_     = calculatedCache_.back();
        toCalculate_    = toCalculateCache_.back();
        toSubtract_     = toSubtractCache_.back();
        removedParents_ = removedParentsCache_.back();
        addedParents_   = addedParentsCache_.back();

        calculatedCache_.pop_back();
        toCalculateCache_.pop_back();
        toSubtractCache_.pop_back();
        removedParentsCache_.pop_back();
        addedParentsCache_.pop_back();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::postAcceptState() {

        //  Clear caches
        calculatedCache_.clear();
        toCalculateCache_.clear();
        toSubtractCache_.clear();
        addedParentsCache_.clear();
        removedParentsCache_.clear();

        addedParents_.clear();
        removedParents_.clear();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::calculateParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps) {
        return ntp_->calculateLogLikelihood(child_, ps);
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::peek() noexcept {

#ifndef NDEBUG
        if (this->value_ <= -std::numeric_limits<Likelihood>::infinity()) {
            std::cerr << "Saved States: " << this->saved_states_stack_.size() << std::endl;
            for (const auto& savedState : this->saved_states_stack_) {
                std::cerr << "State: " << savedState.saved_state << " " << savedState.saved_state_id << std::endl;
            }
        }
#endif
        return Computation::peek();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessV2Impl, typename InfectionEventImpl, typename ParentSetImpl>
    ParentSetDist<InfectionEventImpl> OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessV2Impl, InfectionEventImpl, ParentSetImpl>::calcParentSetDist() {
        ParentSetDist<InfectionEventImpl> dist{};
        auto tmpPs = parentSet_->value();
        tmpPs.insert(latentParent_);
        const int totalNodes = tmpPs.size();
        core::containers::ParentSet<InfectionEventImpl> ps{};
        dist.totalLlik = -std::numeric_limits<Likelihood>::infinity();
        for (int i = 1; i <= ParentSetMaxCardinality; i++) {
            cs_.reset(totalNodes, i);
            while (!cs_.completed) {
                ps.clear();
                for (const auto& idx : cs_.curr) {
                    ps.insert(tmpPs.begin()[idx]);
                }
                Likelihood llik = calculateParentSetLogLikelihood(ps);
                dist.totalLlik  = core::utils::logSumExp(dist.totalLlik, llik);
                dist.parentSetLliks.push_back(std::make_pair(llik, ps));
                cs_.next();
            }
        }

        return dist;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV2_H
