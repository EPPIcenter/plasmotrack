//
// Created by Maxwell Murphy on 2/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESSLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESSLIKELIHOOD_H

#include <boost/math/special_functions/binomial.hpp>
#include <Eigen/Core>

#include "core/datatypes/TransitionMatrix.h"
#include "core/computation/transmission_process/ZTMultiplicativeBinomial.h"



template <int MAX_COI>
class TransmissionProcessLikelihood {
public:

};

//COITransitionMatrix NoSuperInfection::zt_multiplicative_binomial(const NoSuperInfectionParameters &params) {
//    // Calculate the probability of all possible transitions according to a zero-truncated multiplicative binomial
//    // distribution.
//    auto tmp = combinations_matrix.array() *
//               Eigen::pow(params.transmission_complexity_prob_, matrix_one.array()) *
//               Eigen::pow(1 - params.transmission_complexity_prob_, matrix_two.array()) *
//               Eigen::pow(params.transmission_complexity_assoc_, matrix_three.array());
//
//    // Normalize the rows so that it is a valid transition matrix
//    auto res = tmp.array().colwise() / tmp.rowwise().sum();
//    return res;
//}



#endif //TRANSMISSION_NETWORKS_APP_TRANSMISSIONPROCESSLIKELIHOOD_H
