//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include <boost/container/flat_set.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/computation/transmission_process/OrderDerivedParentSet.h"

#include "core/utils/numerics.h"


template <typename TransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename TransmissionEventType>
class OrderBasedTransmissionProcess : public Computation<double>,
                                      public Observable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType>>,
                                      public Cacheable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType>>,
                                      public Checkpointable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType>, double> {
public:
    //// Source Transmission Process changes -> recalculate completely
    //// Node Transmission Process changes -> recalculate completely
    //// Child updated -> recalculate completely
    //// Parent Set adds node -> add node combo to update set, register listener
    //// Parent Set removes node -> subtract node combos immediately, remove listener
    //// Node updated in parent set -> subtract node combos immediately, add node combo to update set
    OrderBasedTransmissionProcess(TransmissionProcessImpl &tp, SourceTransmissionProcessImpl &stp,
            TransmissionEventType &child, OrderDerivedParentSet<TransmissionEventType> &ps) :
            tp_(tp), stp_(stp), ps_(ps), child_(child) {

        tp_.add_set_dirty_listener([&]() {
            transmissionProcessSetDirty();
        });
        tp_.registerCheckpointTarget(*this);

        stp_.add_set_dirty_listener([&]() {
            sourceTransmissionProcessSetDirty();
        });
        stp_.registerCheckpointTarget(*this);

        child_.add_set_dirty_listener([&]() {
            childSetDirty();
        });
        child_.registerCheckpointTarget(*this);

        for(auto& parent : ps_.value()) {


        }

    };



    double value() override {
        return 0;
    }

private:
    friend class Cacheable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType>>;
    friend class Checkpointable<OrderBasedTransmissionProcess<TransmissionProcessImpl, SourceTransmissionProcessImpl, TransmissionEventType>, double>;


    void transmissionProcessSetDirty() {
        this->value_ = 0.0;
        toCalculate = ps_.value();
        this->setDirty();
    };

    void sourceTransmissionProcessSetDirty() {
        this->value_ = 0.0;
        toCalculate = ps_.value();
        this->setDirty();
    };

    void childSetDirty() {
        this->value_ = 0.0;
        toCalculate = ps_.value();
        this->setDirty();
    }


    void removeNode() {};
    void updateNode() {};
    void queueUpdateNode() {};

//    void restoreState() {};

    boost::container::flat_map<TransmissionEventType*, ListenerId_t> preChangeListenerIdMap{};
    boost::container::flat_map<TransmissionEventType*, ListenerId_t> postChangeListenerIdMap{};
    ParentSet<TransmissionEventType> toCalculate{};

    TransmissionProcessImpl &tp_;
    SourceTransmissionProcessImpl &stp_;
    TransmissionEventType &child_;
    OrderDerivedParentSet<TransmissionEventType> &ps_;


};
#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
