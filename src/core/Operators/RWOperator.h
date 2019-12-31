//
// Created by Maxwell Murphy on 12/16/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_RWOPERATOR_H
#define TRANSMISSION_NETWORKS_APP_RWOPERATOR_H

#include "Operator.h"

template <typename T>
struct RWProposal : public Proposal<T> {
    constexpr static float hastings_ratio = 0;
};

template <typename T>
class RWOperator : public Operator<T> {
public:
    RWProposal<T> makeProposal();
};

template<typename T>
RWProposal<T> RWOperator<T>::makeProposal() {
    return RWProposal<T>();
}


#endif //TRANSMISSION_NETWORKS_APP_RWOPERATOR_H
