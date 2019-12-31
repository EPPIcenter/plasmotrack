//
// Created by Maxwell Murphy on 12/11/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_OPERATOR_H
#define TRANSMISSION_NETWORKS_APP_OPERATOR_H

template <typename T>
struct Proposal {
    float hastings_ratio;
    T value;
};

template <typename T>
class Operator {

public:
    virtual Proposal<T> makeProposal() = 0;

};
#endif //TRANSMISSION_NETWORKS_APP_OPERATOR_H

