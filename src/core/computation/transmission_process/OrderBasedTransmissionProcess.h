//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/computation/transmission_process/OrderDerivedParentSet.h"

template <typename TransmissionProcessImpl, typename TransmissionEventType>
class OrderBasedTransmissionProcess : public Computation<double>,
                                      public Observable<OrderBasedTransmissionProcess<TransmissionProcessImpl, TransmissionEventType>>,
                                      public Cacheable<OrderBasedTransmissionProcess<TransmissionProcessImpl, TransmissionEventType>>,
                                      public Checkpointable<OrderBasedTransmissionProcess<TransmissionProcessImpl, TransmissionEventType>, double> {


public:



private:

    TransmissionProcessImpl &tp_;
    OrderDerivedParentSet<TransmissionEventType> &ps_;
    TransmissionEventType &child_;


};
#endif //TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESS_H
