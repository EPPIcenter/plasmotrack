//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MATRIX_H
#define TRANSMISSION_NETWORKS_APP_MATRIX_H

#include <Eigen/Core>

namespace transmission_nets::core::datatypes {

    template<typename T, int DIM>
    using SquareMatrix = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;

    template<int MAX_STATES>
    using LogProbabilityTransitionMatrix = SquareMatrix<float, MAX_STATES>;

    template<int MAX_STATES>
    using TransitionMatrix = SquareMatrix<float, MAX_STATES>;

    template<int MAX_STATES>
    using ProbabilityVector = Eigen::Matrix<float, MAX_STATES, 1>;

    using DynamicProbabilityVector = Eigen::Matrix<float, Eigen::Dynamic, 1>;
    using DynamicArray             = Eigen::Array<float, Eigen::Dynamic, 1>;

    template<int MAX_COI>
    using COITransitionMatrix = LogProbabilityTransitionMatrix<MAX_COI>;

    template<int MAX_TRANSMISSIONS, int MAX_COI>
    using TransmissionProbability = std::array<COITransitionMatrix<MAX_COI>, MAX_TRANSMISSIONS>;

}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_MATRIX_H
