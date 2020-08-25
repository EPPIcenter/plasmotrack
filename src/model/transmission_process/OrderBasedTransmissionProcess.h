//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include <boost/container/flat_set.hpp>
#include <boost/range/adaptors.hpp>

#include "core/computation/PartialLikelihood.h"

#include "core/containers/Infection.h"

#include "OrderDerivedParentSet.h"

#include "core/utils/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"
#include "core/utils/io/serialize.h"


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
class OrderBasedTransmissionProcess : public PartialLikelihood {

    using TransmissionEventType = InfectionEventImpl;
    using ListenerIdMap = boost::container::flat_map<TransmissionEventType *, ListenerId_t>;
public:
    //// Source Transmission Process changes -> subtract source, recalculate
    //// Node Transmission Process changes -> recalculate completely
    //// Child updated -> recalculate completely
    //// Parent Set adds node -> add node combo to update set, register listener
    //// Parent Set removes node -> subtract node combos immediately, remove listener
    //// Node updated in parent set -> subtract node combos immediately, add node combo to update set
    OrderBasedTransmissionProcess(NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp,
                                  TransmissionEventType &child,
                                  OrderDerivedParentSet<TransmissionEventType> &parent_set);


    double value() override;

private:

    void nodeTransmissionProcessSetDirty();

    void sourceTransmissionProcessSetDirty();

    void childSetDirty();

    void addParent(InfectionEventImpl *parent);

    void addParentListeners(InfectionEventImpl *parent);

    void removeParent(InfectionEventImpl *parent);

    void removeParentListeners(InfectionEventImpl *parent);

    void parentUpdated(InfectionEventImpl *parent);

    void customSaveState();

    void customAcceptState();

    void customRestoreState();

    double
    calculateParentLogLikelihoodContribution(InfectionEventImpl *parent, ParentSet<InfectionEventImpl> others);

    ListenerIdMap preChangeListenerIdMap{};
    ListenerIdMap postChangeListenerIdMap{};
    ListenerIdMap saveStateListenerIdMap{};
    ListenerIdMap acceptStateListenerIdMap{};
    ListenerIdMap restoreStateListenerIdMap{};

    // Track parent deltas between save and accept/restore
    boost::container::flat_set<InfectionEventImpl*> addedParents_{};
    boost::container::flat_set<InfectionEventImpl*> removedParents_{};


    boost::container::flat_set<InfectionEventImpl*> toCalculate_{};
    boost::container::flat_set<InfectionEventImpl*> calculated_{};
    boost::container::flat_map<InfectionEventImpl*, double> calculatedParentValues_{};

    boost::container::flat_set<InfectionEventImpl*> toCalculateCache_{};
    boost::container::flat_set<InfectionEventImpl*> calculatedCache_{};
    boost::container::flat_map<InfectionEventImpl*, double> calculatedParentValuesCache_{};


    NodeTransmissionProcessImpl &ntp_;
    SourceTransmissionProcessImpl &stp_;
    InfectionEventImpl &child_;
    OrderDerivedParentSet<InfectionEventImpl> &parentSet_;


    ParentSet<InfectionEventImpl> tmpPs_{};
    std::vector<double> parentLikelihoodContribution_{};
    CombinationIndicesGenerator cs_;

};

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::OrderBasedTransmissionProcess(
        NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp, InfectionEventImpl &child,
        OrderDerivedParentSet<InfectionEventImpl> &parent_set) :
        ntp_(ntp), stp_(stp), child_(child), parentSet_(parent_set) {

    ntp_.add_set_dirty_listener([=, this]() { nodeTransmissionProcessSetDirty(); });
    ntp_.add_save_state_listener([=, this]() { customSaveState(); });
    ntp_.add_restore_state_listener([=, this]() { customRestoreState(); });
    ntp_.add_accept_state_listener([=, this]() { customAcceptState(); });

    stp_.add_set_dirty_listener([=, this]() { sourceTransmissionProcessSetDirty(); });
    stp_.add_save_state_listener([=, this]() { customSaveState(); });
    stp_.add_restore_state_listener([=, this]() { customRestoreState(); });
    stp_.add_accept_state_listener([=, this]() { customAcceptState(); });

    child_.add_post_change_listener([=, this]() { childSetDirty(); });
    child_.add_save_state_listener([=, this]() { customSaveState(); });
    child_.add_restore_state_listener([=, this]() { customRestoreState(); });
    child_.add_accept_state_listener([=, this]() { customAcceptState(); });

    parentSet_.add_element_added_listener([=, this](InfectionEventImpl *parent) { addParent(parent); });
    parentSet_.add_element_removed_listener([=, this](InfectionEventImpl *parent) { removeParent(parent); });
    parentSet_.add_save_state_listener([=, this]() { customSaveState(); });
    parentSet_.add_restore_state_listener([=, this]() { customRestoreState(); });
    parentSet_.add_accept_state_listener([=, this]() { customAcceptState(); });


    for (auto &parent : parentSet_.value()) {
        addParent(parent);
    }

    addedParents_.clear();
    this->setDirty();
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::value() {

    if(this->isDirty()) {
        this->value_ = -std::numeric_limits<double>::infinity();

        for (auto &parent: toCalculate_) {
            calculatedParentValues_.at(parent) = calculateParentLogLikelihoodContribution(parent, calculated_);
            calculated_.insert(parent);
        }
        toCalculate_.clear();

        if(calculatedParentValues_.size() > 0) {
            this->value_ = logSumExp(calculatedParentValues_ | boost::adaptors::map_values);
        }


        this->value_ = logSumExp(this->value_, stp_.value());
        this->setClean();
    }

    assert(this->value_ < std::numeric_limits<double>::infinity());

    return this->value_;
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::nodeTransmissionProcessSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }

    toCalculate_.clear();
    calculated_.clear();
    calculatedParentValues_.clear();

    for (auto &p : parentSet_.value()) {
        calculatedParentValues_.insert_or_assign(p, -std::numeric_limits<double>::infinity());
        toCalculate_.insert(p);
    }

}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::sourceTransmissionProcessSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }
}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::childSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }

    toCalculate_.clear();
    calculated_.clear();
    calculatedParentValues_.clear();

    for (auto &p : parentSet_.value()) {
        calculatedParentValues_.insert_or_assign(p, -std::numeric_limits<double>::infinity());
        toCalculate_.insert(p);
    };
}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::addParent(
        InfectionEventImpl *parent) {
    if (!this->isDirty()) {
        this->setDirty();
    }
    addParentListeners(parent);
    addedParents_.insert(parent);
    toCalculate_.insert(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::addParentListeners(
        InfectionEventImpl *parent) {
//    const auto preChangeListenerId = parent->add_pre_change_listener([=, this]() { parentUpdated(parent); });
    const auto postChangeListenerId = parent->add_post_change_listener([=, this]() { parentUpdated(parent); });
    const auto saveStateListenerId = parent->add_save_state_listener([=, this]() { customSaveState(); });
    const auto acceptStateListenerId = parent->add_accept_state_listener([=, this]() { customAcceptState(); });
    const auto restoreStateListenerId = parent->add_restore_state_listener([=, this]() { customRestoreState(); });

//    preChangeListenerIdMap[parent] = preChangeListenerId;
    postChangeListenerIdMap[parent] = postChangeListenerId;
    saveStateListenerIdMap[parent] = saveStateListenerId;
    acceptStateListenerIdMap[parent] = acceptStateListenerId;
    restoreStateListenerIdMap[parent] = restoreStateListenerId;

    calculatedParentValues_.emplace(parent, -std::numeric_limits<double>::infinity());
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::removeParent(
        InfectionEventImpl *parent) {
    if(!this->isDirty()) {
        this->setDirty();
    }

    calculatedParentValues_.erase(parent);

    if (calculated_.contains(parent)) {
        calculated_.erase(parent);
    } else {
        toCalculate_.erase(parent);
    }

    removeParentListeners(parent);
    removedParents_.insert(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::removeParentListeners(
        InfectionEventImpl *parent) {
//    parent->remove_pre_change_listener(preChangeListenerIdMap[parent]);
    parent->remove_post_change_listener(postChangeListenerIdMap[parent]);
    parent->remove_save_state_listener(saveStateListenerIdMap[parent]);
    parent->remove_accept_state_listener(acceptStateListenerIdMap[parent]);
    parent->remove_restore_state_listener(restoreStateListenerIdMap[parent]);
//    preChangeListenerIdMap.erase(parent);
    postChangeListenerIdMap.erase(parent);
    saveStateListenerIdMap.erase(parent);
    acceptStateListenerIdMap.erase(parent);
    restoreStateListenerIdMap.erase(parent);
}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::parentUpdated(
        InfectionEventImpl *parent) {
    if(!this->isDirty()) {
        this->setDirty();
    }

    if(calculated_.contains(parent)) {
        calculatedParentValues_.at(parent) = -std::numeric_limits<double>::infinity();
        calculated_.erase(parent);
        toCalculate_.insert(parent);
    }
}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::calculateParentLogLikelihoodContribution(
        InfectionEventImpl *parent, ParentSet<InfectionEventImpl> others) {
    const auto otherNodesSize = others.size();
    double max_llik = std::numeric_limits<double>::lowest();
    parentLikelihoodContribution_.clear();
    // Calculate the single parent case
    tmpPs_.clear();
    tmpPs_.insert(parent);
    parentLikelihoodContribution_.push_back(ntp_.calculateLogLikelihood(child_, tmpPs_));
    max_llik = std::max(max_llik, parentLikelihoodContribution_.back());
    for (int i = 1; i < ParentSetMaxCardinality - 1; ++i) {
        cs_.reset(otherNodesSize, i);
        while (!cs_.completed) {
            cs_.next();
            tmpPs_.clear();
            for (const auto &idx : cs_.curr) {
                tmpPs_.insert(others.begin()[idx]);
            }
            tmpPs_.insert(parent);

            parentLikelihoodContribution_.push_back(ntp_.calculateLogLikelihood(child_, tmpPs_));
        }
    }
    double val = logSumExpKnownMax(parentLikelihoodContribution_.begin(), parentLikelihoodContribution_.end(), max_llik);
    assert(val < std::numeric_limits<double>::infinity());
    return val;
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::customSaveState() {
    if (!this->isSaved()) {
        this->notify_save_state();
        this->saved_state_ = this->value();

        calculatedParentValuesCache_ = calculatedParentValues_;
        toCalculateCache_ = toCalculate_;
        calculatedCache_ = calculated_;

    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::customRestoreState() {
    if (this->isSaved()) {
        this->notify_restore_state();
        this->value_ = *(this->saved_state_);

        calculatedParentValues_ = calculatedParentValuesCache_;
        toCalculate_ = toCalculateCache_;
        calculated_ = calculatedCache_;

        for (const auto& parent : addedParents_) {
            removeParentListeners(parent);
        }

        for (const auto& parent : removedParents_) {
            addParentListeners(parent);
        }

        removedParents_.clear();
        addedParents_.clear();

        this->saved_state_.reset();
        this->value();
        this->setClean();
    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::customAcceptState() {
    if (this->isSaved()) {
        this->notify_accept_state();

        removedParents_.clear();
        addedParents_.clear();

        this->saved_state_.reset();
        this->value();
        this->setClean();

    }
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
