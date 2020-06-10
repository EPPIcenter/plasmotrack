//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include <boost/container/flat_set.hpp>
#include <boost/range/adaptors.hpp>

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
    boost::container::flat_set<InfectionEventImpl*> added_parents{};
    boost::container::flat_set<InfectionEventImpl*> removed_parents{};


    boost::container::flat_set<InfectionEventImpl*> to_calculate_{};
    boost::container::flat_set<InfectionEventImpl*> calculated_{};
    boost::container::flat_map<InfectionEventImpl*, double> calculated_parent_values_{};

    boost::container::flat_set<InfectionEventImpl*> to_calculate_cache_{};
    boost::container::flat_set<InfectionEventImpl*> calculated_cache_{};
    boost::container::flat_map<InfectionEventImpl*, double> calculated_parent_values_cache_{};


    NodeTransmissionProcessImpl &ntp_;
    SourceTransmissionProcessImpl &stp_;
    InfectionEventImpl &child_;
    OrderDerivedParentSet<InfectionEventImpl> &parent_set_;


    ParentSet<InfectionEventImpl> tmp_ps_{};
    std::vector<double> parent_likelihood_contribution_{};
    CombinationIndicesGenerator cs_;

};

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::OrderBasedTransmissionProcess(
        NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp, InfectionEventImpl &child,
        OrderDerivedParentSet<InfectionEventImpl> &parent_set) :
        ntp_(ntp), stp_(stp), child_(child), parent_set_(parent_set) {

    ntp_.add_set_dirty_listener([=]() { nodeTransmissionProcessSetDirty(); });
    ntp_.add_save_state_listener([=]() { customSaveState(); });
    ntp_.add_restore_state_listener([=]() { customRestoreState(); });
    ntp_.add_accept_state_listener([=]() { customAcceptState(); });

    stp_.add_set_dirty_listener([=]() { sourceTransmissionProcessSetDirty(); });
    stp_.add_save_state_listener([=]() { customSaveState(); });
    stp_.add_restore_state_listener([=]() { customRestoreState(); });
    stp_.add_accept_state_listener([=]() { customAcceptState(); });

    child_.add_post_change_listener([=]() { childSetDirty(); });
    child_.add_save_state_listener([=]() { customSaveState(); });
    child_.add_restore_state_listener([=]() { customRestoreState(); });
    child_.add_accept_state_listener([=]() { customAcceptState(); });

    parent_set_.add_element_added_listener([=](InfectionEventImpl *parent) { addParent(parent); });
    parent_set_.add_element_removed_listener([=](InfectionEventImpl *parent) { removeParent(parent); });
    parent_set_.add_save_state_listener([=]() { customSaveState(); });
    parent_set_.add_restore_state_listener([=]() { customRestoreState(); });
    parent_set_.add_accept_state_listener([=]() { customAcceptState(); });


    for (auto &parent : parent_set_.value()) {
        addParent(parent);
    }

    added_parents.clear();

    this->setDirty();
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::value() {

    if(this->isDirty()) {
//        this->value_ = -std::numeric_limits<double>::infinity();

        for (auto &parent: to_calculate_) {
            calculated_parent_values_.at(parent) = calculateParentLogLikelihoodContribution(parent, calculated_);
            calculated_.insert(parent);
        }
        to_calculate_.clear();

        if(calculated_parent_values_.size() > 0) {
            this->value_ = logSumExp(calculated_parent_values_ | boost::adaptors::map_values);
        }

//        std::cout << "Total Parents: " << calculated_parent_values_.size() << " -- ";
//        for (const auto& [k, val] : calculated_parent_values_) {
//            std::cout << val << ", ";
//        }

//        std::cout << stp_.value() << std::endl;

        this->value_ = logSumExp(this->value_, stp_.value());
//        this->value_ = logSumExp(this->value_, -1000);
//        std::cout << "Total Value: " << this->value_ << std::endl;
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

    to_calculate_.clear();
    calculated_.clear();
    calculated_parent_values_.clear();

    for (auto &p : parent_set_.value()) {
        calculated_parent_values_.insert_or_assign(p, -std::numeric_limits<double>::infinity());
        to_calculate_.insert(p);
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

    to_calculate_.clear();
    calculated_.clear();
    calculated_parent_values_.clear();

    for (auto &p : parent_set_.value()) {
        calculated_parent_values_.insert_or_assign(p, -std::numeric_limits<double>::infinity());
        to_calculate_.insert(p);
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
    added_parents.insert(parent);
    to_calculate_.insert(parent);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::addParentListeners(
        InfectionEventImpl *parent) {
//    const auto preChangeListenerId = parent->add_pre_change_listener([=]() { parentUpdated(parent); });
    const auto postChangeListenerId = parent->add_post_change_listener([=]() { parentUpdated(parent); });
    const auto saveStateListenerId = parent->add_save_state_listener([=]() { customSaveState(); });
    const auto acceptStateListenerId = parent->add_accept_state_listener([=]() { customAcceptState(); });
    const auto restoreStateListenerId = parent->add_restore_state_listener([=]() { customRestoreState(); });

//    preChangeListenerIdMap[parent] = preChangeListenerId;
    postChangeListenerIdMap[parent] = postChangeListenerId;
    saveStateListenerIdMap[parent] = saveStateListenerId;
    acceptStateListenerIdMap[parent] = acceptStateListenerId;
    restoreStateListenerIdMap[parent] = restoreStateListenerId;

    calculated_parent_values_.emplace(parent, -std::numeric_limits<double>::infinity());
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::removeParent(
        InfectionEventImpl *parent) {
    if(!this->isDirty()) {
        this->setDirty();
    }

    calculated_parent_values_.erase(parent);

    if (calculated_.contains(parent)) {
        calculated_.erase(parent);
    } else {
        to_calculate_.erase(parent);
    }

    removeParentListeners(parent);
    removed_parents.insert(parent);
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
        calculated_parent_values_.at(parent) = -std::numeric_limits<double>::infinity();
        calculated_.erase(parent);
        to_calculate_.insert(parent);
    }
}


template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
double
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::calculateParentLogLikelihoodContribution(
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
    return logSumExpKnownMax(parent_likelihood_contribution_.begin(), parent_likelihood_contribution_.end(), max_llik);
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::customSaveState() {
    if (!this->isSaved()) {
        this->notify_save_state();
        this->saved_state_ = this->value();

        calculated_parent_values_cache_ = calculated_parent_values_;
        to_calculate_cache_ = to_calculate_;
        calculated_cache_ = calculated_;

    }
}

template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
void
OrderBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::customRestoreState() {
    if (this->isSaved()) {
        this->notify_restore_state();
        this->value_ = *(this->saved_state_);

        calculated_parent_values_ = calculated_parent_values_cache_;
        to_calculate_ = to_calculate_cache_;
        calculated_ = calculated_cache_;

        for (const auto& parent : added_parents) {
            removeParentListeners(parent);
        }

        for (const auto& parent : removed_parents) {
            addParentListeners(parent);
        }

        removed_parents.clear();
        added_parents.clear();

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

        removed_parents.clear();
        added_parents.clear();

        this->saved_state_.reset();
        this->value();
        this->setClean();

    }
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
