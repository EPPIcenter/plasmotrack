//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MATRIX_H
#define TRANSMISSION_NETWORKS_APP_MATRIX_H

#include <Eigen/Core>

template <typename T, int DIM>
using SquareMatrix = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;

template <int MAX_STATES>
using ProbabilityMatrix = SquareMatrix<double, MAX_STATES>;


#endif //TRANSMISSION_NETWORKS_APP_MATRIX_H
