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
#include "core/utils/numerics.h"


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

    boost::container::flat_set<InfectionEventImpl*> to_calculate_{};
    boost::container::flat_set<InfectionEventImpl*> calculated_{};

    NodeTransmissionProcessImpl &ntp_;
    SourceTransmissionProcessImpl &stp_;
    InfectionEventImpl &child_;
    OrderDerivedParentSet<InfectionEventImpl> &parent_set_;

    bool source_dirty_ = true;
    double parent_value_cache_ = 0;
    double source_value_cache_ = 0;
    ParentSet<InfectionEventImpl> tmp_ps_{};
    std::vector<double> parent_likelihood_contribution_{};
    CombinationIndicesGenerator cs_;

};

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::OrderBasedTransmissionProcess(
        NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp, InfectionEventImpl &child,
        OrderDerivedParentSet<InfectionEventImpl> &parent_set) :
        ntp_(ntp), stp_(stp), child_(child), parent_set_(parent_set) {

    value_ = 0;
    parent_value_cache_ = 0;
    source_value_cache_ = 0;
    ntp_.add_set_dirty_listener([=]() { nodeTransmissionProcessSetDirty(); });
    ntp_.registerCacheableCheckpointTarget(this);

    stp_.add_set_dirty_listener([=]() { sourceTransmissionProcessSetDirty(); });
    stp_.registerCacheableCheckpointTarget(this);

    child_.add_post_change_listener([=]() { childSetDirty(); });
    child_.registerCacheableCheckpointTarget(this);

    parent_set_.add_element_added_listener([=](InfectionEventImpl *parent) { addParent(parent); });
    parent_set_.add_element_removed_listener([=](InfectionEventImpl *parent) { removeParent(parent); });
    parent_set_.registerCacheableCheckpointTarget(this);

    for (auto &parent : parent_set_.value()) {
        addParent(parent);
        parent->registerCacheableCheckpointTarget(this);
    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::value() {
    if (source_dirty_) {
        source_value_cache_ += exp(stp_.value());
        source_dirty_ = false;
        assert(source_value_cache_ < std::numeric_limits<double>::infinity());
    }

    for (auto &parent: to_calculate_) {
        parent_value_cache_ += calculateParentLikelihoodContribution(parent, calculated_);
        calculated_.insert(parent);
        assert(parent_value_cache_ < std::numeric_limits<double>::infinity());
    }

    to_calculate_.clear();
    this->value_ = parent_value_cache_ + source_value_cache_;
    this->setClean();

    assert(this->value_ < std::numeric_limits<double>::infinity());
//    assert(this->value_ >= 0);

    return this->value_;
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::nodeTransmissionProcessSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }

    this->parent_value_cache_ = 0;

    to_calculate_.clear();
    calculated_.clear();

    for (auto &p : parent_set_.value()) {
        to_calculate_.insert(p);
    };
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::childSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }

    this->parent_value_cache_ = 0;
    this->source_value_cache_ = 0;

    source_dirty_ = true;

    to_calculate_.clear();
    calculated_.clear();

    for (auto &p : parent_set_.value()) {
        to_calculate_.insert(p);
    };
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::sourceTransmissionProcessSetDirty() {
    if(!this->isDirty()) {
        this->setDirty();
    }

    if(!source_dirty_) {
        this->source_value_cache_ = 0;
        source_dirty_ = true;
    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::addParent(
        InfectionEventImpl *parent) {
    if(!this->isDirty()) {
        this->setDirty();
    }
    const auto preChangeListenerId = parent->add_pre_change_listener([=]() { parentUpdated(parent); });
//    const auto postChangeListenerId = parent->add_post_change_listener([=]() { parentUpdated(parent); });
    const auto [saveStateListenerId, acceptStateListenerId, restoreStateListenerId] = parent->registerCacheableCheckpointTarget(this);
    preChangeListenerIdMap[parent] = preChangeListenerId;
//    postChangeListenerIdMap[parent] = postChangeListenerId;
    saveStateListenerIdMap[parent] = saveStateListenerId;
    acceptStateListenerIdMap[parent] = acceptStateListenerId;
    restoreStateListenerIdMap[parent] = restoreStateListenerId;
    to_calculate_.insert(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::removeParent(
        InfectionEventImpl *parent) {
    if(!this->isDirty()) {
        this->setDirty();
    }

    if (calculated_.contains(parent)) {
        calculated_.erase(parent);
        this->parent_value_cache_ -= calculateParentLikelihoodContribution(parent, calculated_);
    } else {
        to_calculate_.erase(parent);
    }

    parent->remove_pre_change_listener(preChangeListenerIdMap[parent]);
//    parent->remove_post_change_listener(postChangeListenerIdMap[parent]);
    parent->remove_save_state_listener(saveStateListenerIdMap[parent]);
    parent->remove_accept_state_listener(acceptStateListenerIdMap[parent]);
    parent->remove_restore_state_listener(restoreStateListenerIdMap[parent]);
    preChangeListenerIdMap.erase(parent);
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
        calculated_.erase(parent);
        parent_value_cache_ -= calculateParentLikelihoodContribution(parent, calculated_);
        to_calculate_.insert(parent);
    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::calculateParentLikelihoodContribution(
        InfectionEventImpl *parent, ParentSet<InfectionEventImpl> others) {
    const auto otherNodesSize = others.size();
    double max_llik = std::numeric_limits<double>::lowest();
    parent_likelihood_contribution_.clear();
    // Calculate the single parent case
    tmp_ps_.clear();
    tmp_ps_.insert(parent);
    parent_likelihood_contribution_.push_back(ntp_.calculateLogLikelihood(child_, tmp_ps_));
    max_llik = std::max(max_llik, parent_likelihood_contribution_.back());
    for (int i = 1; i < ParentSetMaxCardinality - 1; ++i) {
        cs_.reset(otherNodesSize, i);
        while (!cs_.completed) {
            cs_.next();
            tmp_ps_.clear();
            for (const auto &idx : cs_.curr) {
                tmp_ps_.insert(others.begin()[idx]);
            }
            tmp_ps_.insert(parent);

            parent_likelihood_contribution_.push_back(ntp_.calculateLogLikelihood(child_, tmp_ps_));
        }
    }
    return exp(logSumExpKnownMax(parent_likelihood_contribution_.begin(), parent_likelihood_contribution_.end(), max_llik));
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
