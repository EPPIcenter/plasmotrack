//
// Created by Maxwell Murphy on 12/16/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_GIBBSOPERATOR_H
#define TRANSMISSION_NETWORKS_APP_GIBBSOPERATOR_H

#include "Operator.h"

template <typename T>
struct GibbsProposal : public Proposal<T> {
    constexpr static float hastings_ratio = std::numeric_limits::infinity();
};

template <typename T>
class GibbsOperator : public Operator<T> {
    GibbsProposal<T> makeProposal();
};

template<typename T>
GibbsProposal<T> GibbsOperator<T>::makeProposal() {
    return nullptr;
}

#endif //TRANSMISSION_NETWORKS_APP_GIBBSOPERATOR_H
