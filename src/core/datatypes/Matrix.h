//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MATRIX_H
#define TRANSMISSION_NETWORKS_APP_MATRIX_H

#include <Eigen/Core>

namespace transmission_nets::core::datatypes {

    template <typename T, int DIM>
    using SquareMatrix = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;

    template <int MAX_STATES>
    using LogProbabilityTransitionMatrix = SquareMatrix<double, MAX_STATES>;

    template <int MAX_STATES>
    using TransitionMatrix = SquareMatrix<double, MAX_STATES>;

    template <int MAX_STATES>
    using ProbabilityVector = Eigen::Matrix<double, MAX_STATES, 1>;

    using DynamicProbabilityVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using DynamicArray = Eigen::Array<double, Eigen::Dynamic, 1>;

    template <int MAX_COI>
    using COITransitionMatrix = LogProbabilityTransitionMatrix<MAX_COI>;

    template <int MAX_TRANSMISSIONS, int MAX_COI>
    using TransmissionProbability = std::array<COITransitionMatrix<MAX_COI>, MAX_TRANSMISSIONS>;

}


#endif //TRANSMISSION_NETWORKS_APP_MATRIX_H
