//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include <boost/container/flat_set.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/computation/PartialLikelihood.h"

#include "core/containers/Infection.h"

#include "OrderDerivedParentSet.h"

#include "core/utils/CombinationIndicesGenerator.h"


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

    void removeParent(InfectionEventImpl *parent);

    void parentUpdated(InfectionEventImpl *parent);

    double
    calculateParentLikelihoodContribution(InfectionEventImpl *parent, ParentSet<InfectionEventImpl> others);

    ListenerIdMap preChangeListenerIdMap{};
    ListenerIdMap postChangeListenerIdMap{};
    ListenerIdMap saveStateListenerIdMap{};
    ListenerIdMap acceptStateListenerIdMap{};
    ListenerIdMap restoreStateListenerIdMap{};

    std::vector<InfectionEventImpl*> to_calculate_{};
    boost::container::flat_set<InfectionEventImpl*> calculated_{};

    NodeTransmissionProcessImpl &ntp_;
    SourceTransmissionProcessImpl &stp_;
    InfectionEventImpl &child_;
    OrderDerivedParentSet<InfectionEventImpl> &parent_set_;

    bool recalculate_source = true;

};

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::OrderBasedTransmissionProcess(
        NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp, InfectionEventImpl &child,
        OrderDerivedParentSet<InfectionEventImpl> &parent_set) :
        ntp_(ntp), stp_(stp), child_(child), parent_set_(parent_set) {

    value_ = 0;
    ntp_.add_set_dirty_listener([&]() { nodeTransmissionProcessSetDirty(); });
    ntp_.registerCacheableCheckpointTarget(*this);

    stp_.add_set_dirty_listener([&]() { sourceTransmissionProcessSetDirty(); });
    stp_.registerCacheableCheckpointTarget(*this);

    child_.add_post_change_listener([&]() { childSetDirty(); });
    child_.registerCacheableCheckpointTarget(*this);

    parent_set_.add_element_added_listener([&](InfectionEventImpl *parent) { addParent(parent); });
    parent_set_.add_element_removed_listener([&](InfectionEventImpl *parent) { removeParent(parent); });
    parent_set_.registerCacheableCheckpointTarget(*this);

    for (auto &parent : parent_set_.value()) {
        addParent(parent);
    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::value() {
    if (recalculate_source) {
        this->value_ += stp_.value();
        recalculate_source = false;
    }

    for (const auto &parent: to_calculate_) {
        auto parent_likelihood = calculateParentLikelihoodContribution(parent, calculated_);
        this->value_ += parent_likelihood;
        calculated_.insert(parent);
    }
    to_calculate_.clear();
    this->setClean();

    return log(this->value_);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::nodeTransmissionProcessSetDirty() {
    this->setDirty();
    this->value_ = 0;
    recalculate_source = true;
    to_calculate_.clear();
    for (InfectionEventImpl* p : parent_set_.value()) {
        to_calculate_.push_back(p);
    };
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::sourceTransmissionProcessSetDirty() {
    // notify dependencies dirty so may be subtracted, then subtract out after notifying
    this->setDirty();
    this->value_ -= stp_.peek();
    recalculate_source = true;
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::childSetDirty() {
    this->setDirty();
    this->value_ = 0;
    recalculate_source = true;
    to_calculate_.clear();
    for (InfectionEventImpl* p : parent_set_.value()) {
        to_calculate_.push_back(p);
    };
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::addParent(
        InfectionEventImpl *parent) {
//    const auto preChangeListenerId = parent->add_pre_change_listener([&]() { parentUpdated(parent); });
    const auto postChangeListenerId = parent->add_post_change_listener([&]() { parentUpdated(parent); });
    const auto [saveStateListenerId, acceptStateListenerId, restoreStateListenerId] = parent->registerCacheableCheckpointTarget(
            *this);
    this->setDirty();
    postChangeListenerIdMap[parent] = postChangeListenerId;
    saveStateListenerIdMap[parent] = saveStateListenerId;
    acceptStateListenerIdMap[parent] = acceptStateListenerId;
    restoreStateListenerIdMap[parent] = restoreStateListenerId;
    to_calculate_.push_back(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::removeParent(
        InfectionEventImpl *parent) {
    this->setDirty();
    parent->remove_post_change_listener(postChangeListenerIdMap[parent]);
    parent->remove_save_state_listener(saveStateListenerIdMap[parent]);
    parent->remove_accept_state_listener(acceptStateListenerIdMap[parent]);
    parent->remove_restore_state_listener(restoreStateListenerIdMap[parent]);
    calculated_.erase(parent);
    this->value_ -= calculateParentLikelihoodContribution(parent, calculated_);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::parentUpdated(
        InfectionEventImpl *parent) {
    this->setDirty();
    calculated_.erase(parent);
    this->value_ -= calculateParentLikelihoodContribution(parent, calculated_);
    to_calculate_.push_back(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::calculateParentLikelihoodContribution(
        InfectionEventImpl *parent, ParentSet<InfectionEventImpl> others) {
    using ParentSet = ParentSet<InfectionEventImpl>;
    double lik = 0.0;
    const auto otherNodesSize = others.size();
    ParentSet tmp_ps;
    ParentSet p{parent};
    lik += ntp_.calculateLikelihood(child_, p);
    for (int i = 1; i < ParentSetMaxCardinality - 1; ++i) {
        CombinationIndicesGenerator cs(otherNodesSize, i);
        while (!cs.completed) {
            auto indices = cs.next();

            tmp_ps.clear();
            for (const auto &idx : indices) {
                tmp_ps.insert(others.begin()[idx]);
            }
            tmp_ps.insert(parent);

            lik += ntp_.calculateLikelihood(child_, tmp_ps);
        }
    }
    return lik;
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
