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
        OrderBasedTransmissionProcessV2(std::shared_ptr<NodeTransmissionProcessImpl> ntp, std::shared_ptr<SourceTransmissionProcessImpl> stp, std::shared_ptr<InfectionEventImpl> child, std::shared_ptr<ParentSetImpl> parent_set, std::shared_ptr<InfectionEventImpl> latent_parent);

        Likelihood value() override;

        Likelihood recalculate(bool verbose = false);

        void initializeValue();

        Likelihood calculateParentLogLikelihoodContribution(std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others);
        Likelihood peekParentLogLikelihoodContribution(std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others);

        Likelihood calculateLatentParentLogLikelihoodContribution(const core::containers::ParentSet<InfectionEventImpl>& others);
        Likelihood peekLatentParentLogLikelihoodContribution(const core::containers::ParentSet<InfectionEventImpl>& others);

        Likelihood calculateParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps);
        Likelihood peekParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps);

        ParentSetDist<InfectionEventImpl> calcParentSetDist();

        std::string identifier() override;

        std::shared_ptr<NodeTransmissionProcessImpl> ntp_;
        std::shared_ptr<SourceTransmissionProcessImpl> stp_;
        std::shared_ptr<InfectionEventImpl> child_;
        std::shared_ptr<ParentSetImpl> parentSet_;
        std::shared_ptr<InfectionEventImpl> latentParent_;

    private:
        void nodeTransmissionProcessSetDirty();

        void sourceTransmissionProcessSetDirty();

        void childSetDirty();

        void addParent(std::shared_ptr<InfectionEventImpl> parent);

        void addParentListeners(std::shared_ptr<InfectionEventImpl> parent);

        void addLatentParentListeners();

        void removeParent(std::shared_ptr<InfectionEventImpl> parent);

        void removeParentListeners(std::shared_ptr<InfectionEventImpl> parent);

        void preParentUpdated(std::shared_ptr<InfectionEventImpl> parent);

        void postParentUpdated(std::shared_ptr<InfectionEventImpl> parent);

        void preLatentParentUpdated();

        void postLatentParentUpdated();

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
        bool latentParentDirty_ = true;

        std::deque<InfectionEventSet> calculatedCache_{};
        std::deque<InfectionEventSet> toCalculateCache_{};
        std::deque<std::vector<Likelihood>> toSubtractCache_{};

        // helper vars for calculating parent likelihood contributions
        core::containers::ParentSet<InfectionEventImpl> tmpPs_{};
        std::vector<Likelihood> parentLikelihoodContribution_{};
        std::vector<Likelihood> toAdd_{};

        std::string action_;

    public:
        Likelihood peek() noexcept override;
    };


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::OrderBasedTransmissionProcessV2(
            std::shared_ptr<NodeTransmissionProcessImpl> ntp, std::shared_ptr<SourceTransmissionProcessImpl> stp, std::shared_ptr<InfectionEventImpl> child, std::shared_ptr<ParentSetImpl> parent_set, std::shared_ptr<InfectionEventImpl> latent_parent) : ntp_(std::move(ntp)), stp_(std::move(stp)), child_(std::move(child)), parentSet_(std::move(parent_set)), latentParent_(std::move(latent_parent)) {

        ntp_->add_set_dirty_listener([=, this]() { nodeTransmissionProcessSetDirty(); });
        ntp_->registerCacheableCheckpointTarget(this);

        stp_->add_set_dirty_listener([=, this]() { sourceTransmissionProcessSetDirty(); });
        stp_->registerCacheableCheckpointTarget(this);

        child_->add_post_change_listener([=, this]() { childSetDirty(); });
        child_->registerCacheableCheckpointTarget(this);

        parentSet_->add_element_added_listener([=, this](std::shared_ptr<InfectionEventImpl> parent) { addParent(parent); });
        parentSet_->add_element_removed_listener([=, this](std::shared_ptr<InfectionEventImpl> parent) { removeParent(parent); });
        parentSet_->registerCacheableCheckpointTarget(this);

        this->addPostSaveHook([=, this]() { this->postSaveState(); });
        this->addPostRestoreHook([=, this]() { this->postRestoreState(); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        addLatentParentListeners();

        // Register the latent parent and the parents in the derived parent set
        for (auto& parent : parentSet_->value()) {
            addParent(parent);
        }

        // Clear the added parents tracker on first initialization
        addedParents_.clear();
        this->value_ = -std::numeric_limits<Likelihood>::infinity();
        this->setDirty();
        this->initializeValue();
//        this->value();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    std::string OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::identifier() {
        return "OrderBasedTransmissionProcessV2";
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::value() {
        if (this->isDirty()) {
            [[maybe_unused]] auto preValue = this->value_;
            toAdd_.clear();
//            Likelihood maxLlik = -std::numeric_limits<Likelihood>::infinity();
//            Likelihood totalToSubtract = 0.0;
//            if (!toSubtract_.empty()) {
//                totalToSubtract = core::utils::logSumExp(toSubtract_);
//            }
            this->value_ = recalculate();

//            if (std::abs(this->value_ - totalToSubtract) < 20) {
//                // numerical precision issues when subtracting two very close numbers
//                this->value_ = recalculate();
//            } else {
//                if (totalToSubtract > 1e-6) {
//                    this->value_ = core::utils::logDiffExp(this->value_, totalToSubtract);
//                }
//
//                toAdd_.push_back(this->value_);
//                maxLlik = std::max(maxLlik, toAdd_.back());
//
//                // If the latent parent is dirty, recalculate all parent sets including the latent parent
//                if (this->latentParentDirty_) {
//                    if (calculated_.size() > 0) {
//                        toAdd_.push_back(calculateLatentParentLogLikelihoodContribution(calculated_));
//                        maxLlik = std::max(maxLlik, toAdd_.back());
//                    }
//                    toAdd_.push_back(ntp_->calculateLogLikelihood(child_, latentParent_, stp_));
//                    maxLlik = std::max(maxLlik, toAdd_.back());
//                }
//
//                // then recalculate all parent sets that are dirty (automatically includes the latent parent)
//                for (const auto& parent : toCalculate_) {
//                    assert(!(calculated_.contains(parent)));
//                    // calculate the needed parent sets, including the latent parent
//                    toAdd_.push_back(calculateParentLogLikelihoodContribution(parent, calculated_));
//                    maxLlik = std::max(maxLlik, toAdd_.back());
//                    calculated_.insert(parent);
//                }
//
//
//                this->value_ = core::utils::logSumExpKnownMax(toAdd_.begin(), toAdd_.end(), maxLlik);
//
//                assert(this->value_ < std::numeric_limits<Likelihood>::infinity());
//
//#ifdef DEBUG_LIKELIHOOD
////
////                auto tmp = recalculate();
////                if (std::abs(tmp - this->value_) > .1) {
////                    fmt::print("Action: {}", action_);
////                    fmt::print("tmp value: {}, {}, {}\n", tmp, preValue, this->value_);
//////                if (this->value_ > -std::numeric_limits<Likelihood>::infinity() and (std::abs(tmp - this->value_) > 1)) {
//////                    fmt::print("Err: {}, {}, {}, {} \n", std::abs(tmp - this->value_), tmp, this->value_, preValue);
//////                    fmt::print("{} Likelihood mismatch OBTP: {}, {}\n", fmt::ptr(this), tmp, this->value_);
////                    auto second_validate = recalculate(true);
////                    fmt::print("Second Validation: {}\n", second_validate);
////                    fmt::print("Child Node: {}\n", child_->id());
////                    fmt::print("Total Parents: {}\n", parentSet_->value().size());
////                    fmt::print("To Add: {}\n", fmt::join(toAdd_, ", "));
////                    fmt::print("To Subtract: {}\n", fmt::join(toSubtract_, ", "));
////                    fmt::print("Total to subtract: {0:.36f}\n", totalToSubtract);
////                    std::vector<std::string> removedParentIDs{};
////                    std::vector<std::string> currentParentIDs{};
////                    std::vector<std::string> toCalculateIDs{};
////                    std::vector<std::string> calculatedParents{};
////                    std::transform(removedParents_.begin(), removedParents_.end(), std::back_inserter(removedParentIDs), [](auto p) { return p->id(); });
////                    std::transform(toCalculate_.begin(), toCalculate_.end(), std::back_inserter(toCalculateIDs), [](auto p) { return p->id(); });
////                    std::transform(
////                            calculated_.begin(), calculated_.end(), std::back_inserter(calculatedParents), [](auto p) { return p->id(); });
////                    for (const auto& parent : parentSet_->value()) {
////                        currentParentIDs.push_back(parent->id());
////                    }
////                    fmt::print("Removed Parents: {}\n", fmt::join(removedParentIDs, ", "));
////                    fmt::print("Current Parents: {}\n", fmt::join(currentParentIDs, ", "));
////                    fmt::print("To Calculate: {}\n", fmt::join(toCalculateIDs, ", "));
////                    fmt::print("Calculated Parents: {}\n", fmt::join(calculatedParents, ", "));
////
////                    exit(1);
////                }
//#endif
//            }
//            fmt::print("Pre Value: {}\nPost Value: {}\n", preValue, this->value_);
            this->setClean();
            this->latentParentDirty_ = false;
            toSubtract_.clear();
            toCalculate_.clear();
        }
        return this->value_;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::recalculate(bool verbose) {
        /*
         * function that completely recalculates the likelihood.
         */
        auto val = -std::numeric_limits<Likelihood>::infinity();
        InfectionEventSet tmpCalculated_{};
        std::vector<Likelihood> tmpToAdd_{};
        Likelihood maxLlik = val;

        std::vector<std::string> parentIDs{};

        for (const auto& parent : parentSet_->value()) {
            tmpToAdd_.push_back(calculateParentLogLikelihoodContribution(parent, tmpCalculated_));
            maxLlik = std::max(maxLlik, tmpToAdd_.back());
            tmpCalculated_.insert(parent);
            parentIDs.push_back(parent->id());
        }

        // calculate the single parent set consisting of the latent parent alone
        tmpToAdd_.push_back(ntp_->calculateLogLikelihood(child_, latentParent_, stp_));
        maxLlik = std::max(maxLlik, tmpToAdd_.back());

        val         = core::utils::logSumExpKnownMax(tmpToAdd_.begin(), tmpToAdd_.end(), maxLlik);
        calculated_ = std::move(tmpCalculated_);


//        if (verbose or (parentSet_->value().size() == 0 and val <= -std::numeric_limits<Likelihood>::infinity())) {
        if (verbose) {
            fmt::print("Action: {}\n", action_);
            fmt::print("Child Node: {}\n", child_->id());
            fmt::print("Validate to Add: {}\n", fmt::join(tmpToAdd_, ", "));
            fmt::print("Validation Current Parents: {}\n", fmt::join(parentIDs, ", "));
            fmt::print("Latent parent contribution {}\n", tmpToAdd_.back());
            fmt::print("Validation MaxLlik: {}\n", maxLlik);
            fmt::print("Validation Value: {}\n", val);
            fmt::print("====================================\n");
        }
        return val;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::initializeValue() {
        this->value_ = recalculate(true);
        this->latentParentDirty_ = false;
        this->setClean();
        toSubtract_.clear();
        toCalculate_.clear();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::nodeTransmissionProcessSetDirty() {
        /*
         * Completely recalculate the likelihood
         */

        action_ = "Node Transmission Process Set Dirty";

        if (!this->isDirty()) {
            this->setDirty();
        }

        // recalculate the latentParent parent set
        this->latentParentDirty_ = true;

        toCalculate_.clear();
        toSubtract_.clear();
        calculated_.clear();

        for ([[maybe_unused]] auto& p : parentSet_->value()) {
            toCalculate_.insert(p);
        }

        this->value_ = -std::numeric_limits<Likelihood>::infinity();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::sourceTransmissionProcessSetDirty() {
        /*
         * Recalculate latent parent contribution
         */

        action_ = "Source Transmission Process Set Dirty";

        if (!this->isDirty()) {
            this->setDirty();
        }

        // recalculate the latentParent parent set
        if (!latentParentDirty_) {
            if (calculated_.size() > 0) {
                // Calculate the parent sets that contain the latent parent
                Likelihood latentParentContr = peekLatentParentLogLikelihoodContribution(calculated_);
                toSubtract_.push_back(latentParentContr);
            }
            Likelihood latentParentSetContr = ntp_->peekLogLikelihood(child_, latentParent_, stp_);
            toSubtract_.push_back(latentParentSetContr);
            this->latentParentDirty_ = true;
        }

    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::childSetDirty() {
        /*
         * Completely recalculate the likelihood
         */

        action_ = "Child set dirty";

        if (!this->isDirty()) {
            this->setDirty();
        }

        // recalculate the latentParent parent set
        this->latentParentDirty_ = true;

        toCalculate_.clear();
        toSubtract_.clear();
        calculated_.clear();

        for ([[maybe_unused]] auto& p : parentSet_->value()) {
            toCalculate_.insert(p);
        };
        this->value_ = -std::numeric_limits<Likelihood>::infinity();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::addParent(
            std::shared_ptr<InfectionEventImpl> parent) {
        /*
         * Mark new parent to be calculated, add listeners to parent
         */

        action_ = "Add Parent";

        addParentListeners(parent);
        addedParents_.insert(parent);
        toCalculate_.insert(parent);

        if (!this->isDirty()) {
            this->setDirty();
        }
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::addLatentParentListeners() {
        const auto preChangeListenerId    = latentParent_->add_pre_change_listener([=, this]() { preLatentParentUpdated(); });
        const auto postChangeListenerId   = latentParent_->add_post_change_listener([=, this]() { postLatentParentUpdated(); });
        const auto saveStateListenerId    = latentParent_->add_save_state_listener([=, this](const std::string& savedStateId) { saveState(savedStateId); });
        const auto acceptStateListenerId  = latentParent_->add_accept_state_listener([=, this]() { acceptState(); });
        const auto restoreStateListenerId = latentParent_->add_restore_state_listener([=, this](const std::string& savedStateId) { restoreState(savedStateId); });

        assert(!preChangeListenerIdMap.contains(latentParent_));
        preChangeListenerIdMap[latentParent_]    = preChangeListenerId;
        postChangeListenerIdMap[latentParent_]   = postChangeListenerId;
        saveStateListenerIdMap[latentParent_]    = saveStateListenerId;
        acceptStateListenerIdMap[latentParent_]  = acceptStateListenerId;
        restoreStateListenerIdMap[latentParent_] = restoreStateListenerId;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::addParentListeners(
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


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::removeParent(
            std::shared_ptr<InfectionEventImpl> parent) {

        action_ = "Remove Parent";

        calculated_.erase(parent);

        toSubtract_.push_back(calculateParentLogLikelihoodContribution(parent, calculated_));

        removeParentListeners(parent);
        removedParents_.insert(parent);

        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::removeParentListeners(
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


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::preParentUpdated(
            std::shared_ptr<InfectionEventImpl> parent) {
        // Before the parent is updated, we need to calculate the contribution
        // of the parent from the likelihood, including sets with the latent parent.
        if (calculated_.contains(parent)) {
            calculated_.erase(parent);
            toSubtract_.push_back(peekParentLogLikelihoodContribution(parent, calculated_));
            toCalculate_.insert(parent);
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postParentUpdated(
            [[maybe_unused]] std::shared_ptr<InfectionEventImpl> parent) {
        // After the parent is updated, set the likelihood to dirty.
        action_ = "Parent Updated";
        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::preLatentParentUpdated() {
        // Before the latent parent is updated, we need to calculate the contribution
        // of the latent parent from the likelihood.
        if (!latentParentDirty_) {
            if (calculated_.size() > 0) {
                // Calculate the parent sets that contain the latent parent
                Likelihood latentParentContr = peekLatentParentLogLikelihoodContribution(calculated_);
                toSubtract_.push_back(latentParentContr);
            }
            Likelihood latentParentSetContr = ntp_->peekLogLikelihood(child_, latentParent_, stp_);
            toSubtract_.push_back(latentParentSetContr);
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postLatentParentUpdated() {
        // After the latent parent is updated, set the likelihood and latentParent flag to dirty.
        action_ = "Latent Parent Updated";
        latentParentDirty_ = true;
        if (!this->isDirty()) {
            this->setDirty();
        }
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::calculateParentLogLikelihoodContribution(
            std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others) {
        // Calculate the contribution of the parent to the likelihood. Includes with and without the latent parent.
        core::utils::generators::CombinationIndicesGenerator comboGen;
        const int otherNodesSize = others.size();
        parentLikelihoodContribution_.clear();

        // Calculate the single parent case
        tmpPs_.clear();
        tmpPs_.insert(parent);

        parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, tmpPs_));
        Likelihood max_llik = parentLikelihoodContribution_.back();
        parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, latentParent_, tmpPs_, stp_));
        max_llik = std::max(max_llik, parentLikelihoodContribution_.back());

        for (int i = 1; i < ParentSetMaxCardinality and i <= otherNodesSize; ++i) {
            comboGen.reset(otherNodesSize, i);
            while (!comboGen.completed) {
                // Generate the parent set
                tmpPs_.clear();
                for (const auto& idx : comboGen.curr) {
                    tmpPs_.insert(others.begin()[idx]);
                }
                tmpPs_.insert(parent);

                // calculate without the latent parent
                parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, tmpPs_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());

                // calculate with the latent parent
                parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, latentParent_, tmpPs_, stp_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                comboGen.next();
            }
        }
        Likelihood val = core::utils::logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
        assert(val < std::numeric_limits<Likelihood>::infinity());
        return val;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::peekParentLogLikelihoodContribution(
            std::shared_ptr<InfectionEventImpl> parent, const core::containers::ParentSet<InfectionEventImpl>& others) {

        core::utils::generators::CombinationIndicesGenerator comboGen;
        const int otherNodesSize = others.size();
        parentLikelihoodContribution_.clear();

        // Calculate the single parent case
        tmpPs_.clear();
        tmpPs_.insert(parent);

        parentLikelihoodContribution_.push_back(ntp_->peekLogLikelihood(child_, tmpPs_));
        Likelihood max_llik = parentLikelihoodContribution_.back();
        parentLikelihoodContribution_.push_back(ntp_->peekLogLikelihood(child_, latentParent_, tmpPs_, stp_));
        max_llik = std::max(max_llik, parentLikelihoodContribution_.back());

        for (int i = 1; i < ParentSetMaxCardinality and i <= otherNodesSize; ++i) {
            comboGen.reset(otherNodesSize, i);
            while (!comboGen.completed) {
                // Generate the parent set
                tmpPs_.clear();
                for (const auto& idx : comboGen.curr) {
                    tmpPs_.insert(others.begin()[idx]);
                }
                tmpPs_.insert(parent);

                parentLikelihoodContribution_.push_back(ntp_->peekLogLikelihood(child_, tmpPs_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                parentLikelihoodContribution_.push_back(ntp_->peekLogLikelihood(child_, latentParent_, tmpPs_, stp_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                comboGen.next();
            }
        }
        Likelihood val = core::utils::logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
        assert(val < std::numeric_limits<Likelihood>::infinity());
        return val;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::calculateLatentParentLogLikelihoodContribution(const core::containers::ParentSet<InfectionEventImpl>& others) {
        // Calculate all parent subsets from others including the latent parent.
        core::utils::generators::CombinationIndicesGenerator comboGen;
        const int otherNodesSize = others.size();
        parentLikelihoodContribution_.clear();
        Likelihood max_llik = -std::numeric_limits<Likelihood>::infinity();

        for (int i = 1; i < ParentSetMaxCardinality and i <= otherNodesSize; ++i) {
            comboGen.reset(otherNodesSize, i);
            while (!comboGen.completed) {
                tmpPs_.clear();
                for (const auto& idx : comboGen.curr) {
                    tmpPs_.insert(others.begin()[idx]);
                }
                parentLikelihoodContribution_.push_back(ntp_->calculateLogLikelihood(child_, latentParent_, tmpPs_, stp_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                comboGen.next();
            }
        }

        Likelihood val = core::utils::logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
        assert(val < std::numeric_limits<Likelihood>::infinity());
        return val;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::peekLatentParentLogLikelihoodContribution(const core::containers::ParentSet<InfectionEventImpl>& others) {
        // Calculate all parent subsets from others including the latent parent.
        core::utils::generators::CombinationIndicesGenerator comboGen;
        const int otherNodesSize = others.size();
        parentLikelihoodContribution_.clear();
        Likelihood max_llik = -std::numeric_limits<Likelihood>::infinity();

        for (int i = 1; i < ParentSetMaxCardinality and i <= otherNodesSize; ++i) {
            comboGen.reset(otherNodesSize, i);
            while (!comboGen.completed) {
                tmpPs_.clear();
                for (const auto& idx : comboGen.curr) {
                    tmpPs_.insert(others.begin()[idx]);
                }
                parentLikelihoodContribution_.push_back(ntp_->peekLogLikelihood(child_, latentParent_, tmpPs_, stp_));
                max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
                comboGen.next();
            }
        }

        Likelihood val = core::utils::logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
        assert(val < std::numeric_limits<Likelihood>::infinity());
        return val;
    }




    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postSaveState() {
        calculatedCache_.emplace_back(calculated_);
        toCalculateCache_.emplace_back(toCalculate_);
        toSubtractCache_.emplace_back(toSubtract_);

        addedParentsCache_.emplace_back(addedParents_);
        removedParentsCache_.emplace_back(removedParents_);

        addedParents_.clear();
        removedParents_.clear();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postRestoreState() {
        for ([[maybe_unused]] const auto& parent : addedParents_) {
            removeParentListeners(parent);
        }

        for ([[maybe_unused]] const auto& parent : removedParents_) {
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
        this->setDirty();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postAcceptState() {

        //  Clear caches
        calculatedCache_.clear();
        toCalculateCache_.clear();
        toSubtractCache_.clear();
        addedParentsCache_.clear();
        removedParentsCache_.clear();

        addedParents_.clear();
        removedParents_.clear();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::calculateParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps) {
        return ntp_->calculateLogLikelihood(child_, ps);
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::peekParentSetLogLikelihood(const core::containers::ParentSet<InfectionEventImpl>& ps) {
        return ntp_->peekLogLikelihood(child_, ps);
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::peek() noexcept {

#ifndef NDEBUG
        if (this->value_ <= -std::numeric_limits<Likelihood>::infinity()) {
            std::cerr << "Saved States: " << this->saved_states_stack_.size() << std::endl;
            for ([[maybe_unused]] const auto& savedState : this->saved_states_stack_) {
                std::cerr << "State: " << savedState.saved_state << " " << savedState.saved_state_id << std::endl;
            }
        }
#endif
        return Computation::peek();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    ParentSetDist<InfectionEventImpl> OrderBasedTransmissionProcessV2<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::calcParentSetDist() {
        ParentSetDist<InfectionEventImpl> dist{};
        auto tmpPs = parentSet_->value();
        const int totalNodes = tmpPs.size();
        core::containers::ParentSet<InfectionEventImpl> ps{};
        Likelihood llik = ntp_->calculateLogLikelihood(child_, latentParent_, stp_);
        dist.totalLlik = llik;
        ps.insert(latentParent_);
        dist.parentSetLliks.push_back(std::make_pair(llik, ps));

        core::utils::generators::CombinationIndicesGenerator comboGen;
        for (int i = 1; i <= ParentSetMaxCardinality; i++) {
            comboGen.reset(totalNodes, i);
            while (!comboGen.completed) {
                ps.clear();
                for (const auto& idx : comboGen.curr) {
                    ps.insert(tmpPs.begin()[idx]);
                }

                // calculate parentset without latent parent
                llik = ntp_->calculateLogLikelihood(child_, ps);
                dist.totalLlik  = core::utils::logSumExp(dist.totalLlik, llik);
                dist.parentSetLliks.push_back(std::make_pair(llik, ps));

                // calculate parentset with latent parent
                llik = ntp_->calculateLogLikelihood(child_, latentParent_, ps, stp_);
                ps.insert(latentParent_);
                dist.totalLlik  = core::utils::logSumExp(dist.totalLlik, llik);
                dist.parentSetLliks.push_back(std::make_pair(llik, ps));
                comboGen.next();
            }
        }

        return dist;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV2_H
