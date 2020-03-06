//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include <boost/container/flat_set.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "OrderDerivedParentSet.h"

#include "core/utils/numerics.h"
#include "core/utils/CombinationIndicesGenerator.h"


template<typename TransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename TransmissionEventType, int ParentSetMaxCardinality>
class OrderBasedTransmissionProcess : public Computation<double>,
                                      public Observable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType, ParentSetMaxCardinality>>,
                                      public Cacheable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType, ParentSetMaxCardinality>>,
                                      public Checkpointable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType, ParentSetMaxCardinality>, double> {

    using ListenerIdMap = boost::container::flat_map<TransmissionEventType *, ListenerId_t>;
public:
    //// Source Transmission Process changes -> subtract source, recalculate
    //// Node Transmission Process changes -> recalculate completely
    //// Child updated -> recalculate completely
    //// Parent Set adds node -> add node combo to update set, register listener
    //// Parent Set removes node -> subtract node combos immediately, remove listener
    //// Node updated in parent set -> subtract node combos immediately, add node combo to update set
    OrderBasedTransmissionProcess(TransmissionProcessImpl &tp, SourceTransmissionProcessImpl &stp,
                                  TransmissionEventType &child, OrderDerivedParentSet<TransmissionEventType> &ps) :
            tp_(tp), stp_(stp), ps_(ps), child_(child) {

        value_ = 0;
        tp_.add_set_dirty_listener([&]() { transmissionProcessSetDirty(); });
        tp_.registerCheckpointTarget(*this);

        stp_.add_set_dirty_listener([&]() { sourceTransmissionProcessSetDirty(); });
        stp_.registerCheckpointTarget(*this);

        child_.add_post_change_listener([&]() { childSetDirty(); });
        child_.registerCheckpointTarget(*this);

        ps_.add_element_added_listener([&](TransmissionEventType *parent) { addParent(parent); });
        ps_.add_element_removed_listener([&](TransmissionEventType *parent) { removeParent(parent); });
        ps_.registerCheckpointTarget(*this);

        for (auto &parent : ps_.value()) {
            addParent(parent);
        }

    };


    double value() override {
        if(recalculateSource) {
            this->value_ += stp_.calculateLikelihood(child_);
        }

        for(const auto& parent: toCalculate) {
            this->value_ += calculateParentLikelihoodContribution(parent, calculated);
            calculated.insert(parent);
        }
        toCalculate.clear();

        return this->value_;
    }

private:
    friend class Cacheable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType, ParentSetMaxCardinality>>;

    friend class Checkpointable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType, ParentSetMaxCardinality>, double>;

    void transmissionProcessSetDirty() {
        recalculateSource = true;
        toCalculate = std::vector(ps_.value().begin(), ps_.value().end());
        this->setDirty();
    };

    void sourceTransmissionProcessSetDirty() {
        // notify dependencies dirty so may be subtracted, then subtract out after notifying
        this->setDirty();
        this->value_ -= stp_.peekCalculateLikelihood(child_);
        recalculateSource = true;
        this->setClean();
    };

    void childSetDirty() {
        this->setDirty();
        this->value_ = 0;
        recalculateSource = true;
        toCalculate = std::vector(ps_.value().begin(), ps_.value().end());
    }


    void addParent(TransmissionEventType *parent) {
        const auto preChangeListenerId = parent->add_pre_change_listener([&]() { parentUpdated(parent); });
        const auto postChangeListenerId = parent->add_post_change_listener([&]() { postParentUpdated(parent); });
        const auto[saveStateListenerId, acceptStateListenerId, restoreStateListenerId] = parent->registerCheckpointTarget(
                *this);
        preChangeListenerIdMap[parent] = preChangeListenerId;
        postChangeListenerIdMap[parent] = postChangeListenerId;
        saveStateListenerIdMap[parent] = saveStateListenerId;
        acceptStateListenerIdMap[parent] = acceptStateListenerId;
        restoreStateListenerIdMap[parent] = restoreStateListenerId;
        toCalculate.push_back(parent);
    };

    void removeParent(TransmissionEventType *parent) {
        parent->remove_pre_change_listener(preChangeListenerIdMap[parent]);
        parent->remove_post_change_listener(postChangeListenerIdMap[parent]);
        parent->remove_save_state_listener(saveStateListenerIdMap[parent]);
        parent->remove_accept_state_listener(acceptStateListenerIdMap[parent]);
        parent->remove_restore_state_listener(restoreStateListenerIdMap[parent]);

        calculated.erase(parent);
        this->value_ -= calculateParentLikelihoodContribution(parent, calculated);
    };

    void parentUpdated(TransmissionEventType *parent) {
        this->setDirty();
        calculated.erase(parent);
        this->value_ -= calculateParentLikelihoodContribution(parent, calculated);
        toCalculate.push_back(parent);
    };

    double calculateParentLikelihoodContribution(TransmissionEventType* parent, ParentSet<TransmissionEventType> others) {
        double lik = 0.0;
        const auto otherNodesSize = others.size();
        ParentSet<TransmissionEventType> tmp_ps;
        lik += tp_.calculateLikelihood(child_, ParentSet(parent));
        for (int i = 1; i < ParentSetMaxCardinality - 1; ++i) {
            CombinationIndicesGenerator cs(otherNodesSize, i);
            while (!cs.completed) {
                auto indices = cs.next();

                tmp_ps.clear();
                for(const auto& idx : indices) {
                    tmp_ps.insert(others.begin()[idx]);
                }
                tmp_ps.insert(parent);

                lik += tp_.calculateLikelihood(child_, tmp_ps);
            }
        }
        return lik;
    }

    ListenerIdMap preChangeListenerIdMap{};
    ListenerIdMap postChangeListenerIdMap{};
    ListenerIdMap saveStateListenerIdMap{};
    ListenerIdMap acceptStateListenerIdMap{};
    ListenerIdMap restoreStateListenerIdMap{};

    std::vector<TransmissionEventType> toCalculate{};
    boost::container::flat_set<TransmissionEventType> calculated{};

    TransmissionProcessImpl &tp_;
    SourceTransmissionProcessImpl &stp_;
    TransmissionEventType &child_;
    OrderDerivedParentSet<TransmissionEventType> &ps_;

    bool recalculateSource = false;

};

#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
