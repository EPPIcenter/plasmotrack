//
// Created by Maxwell Murphy on 10/16/19.
//

#ifndef TRANSMISSION_NETWORK_MODEL_TRANSMISSIONMECHANISM_H
#define TRANSMISSION_NETWORK_MODEL_TRANSMISSIONMECHANISM_H

#include <array>
#include <eigen3/Eigen/Core>
#include <boost/math/special_functions/binomial.hpp>

#include "TransmissionNetworkModelConfig.h"
#include "utils.h"
#include "Edge.h"

typedef Eigen::Matrix<Probability, MAX_COI, MAX_COI, Eigen::RowMajor> COITransitionMatrix;
typedef std::array<COITransitionMatrix, MAX_TRANSMISSIONS> TransmissionProbability;
typedef Eigen::Matrix<float, MAX_COI, MAX_COI, Eigen::RowMajor> TransitionMatrix;

template<class TransmissionMechanismParameters>
class TransmissionMechanism {
private:
    TransmissionMechanismParameters transmission_mech_parameters;
public:
//    TransmissionMechanism(const TransmissionMechanismParameters &transmission_mech_parameters);
//    TransmissionMechanism();

//    virtual LogLik calc_llik_transmission(
//            const ObservedNode &parent_node, const ObservedNode &child_node, int num_transmissions,
//            const std::vector<int> &total_alleles, int total_loci
//            ) const=0;
//
//    virtual LogLik calc_llik_transmission(
//            const SourceNode &parent_node, const ObservedNode &child_node, const std::vector<int> &total_alleles, int total_loci
//            ) const=0;

    virtual LogLik calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set, int total_loci,
                                          const std::vector<int> &total_alleles) = 0;

    virtual LogLik calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set,
                                          const SourceNode &source_node, int total_loci,
                                          const std::vector<int> &total_alleles) = 0;

    virtual LogLik calc_llik_transmission(const ObservedNode &node, const SourceNode &source_node, int total_loci,
                                          const std::vector<int> &total_alleles) = 0;

    virtual void set_transmission_mech_parameters(const TransmissionMechanismParameters &params) = 0;

    virtual TransmissionMechanismParameters get_transmission_mech_parameters() = 0;

};



#endif //TRANSMISSION_NETWORK_MODEL_TRANSMISSIONMECHANISM_H
