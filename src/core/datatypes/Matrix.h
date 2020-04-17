//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MATRIX_H
#define TRANSMISSION_NETWORKS_APP_MATRIX_H

#include <Eigen/Core>

template <typename T, int DIM>
using SquareMatrix = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;

template <int MAX_STATES>
using LogProbabilityMatrix = SquareMatrix<double, MAX_STATES>;

template <int MAX_STATES>
using ProbabilityMatrix = SquareMatrix<double, MAX_STATES>;

template <int MAX_STATES>
using ProbabilityVector = Eigen::Matrix<double, MAX_STATES, 1>;

using DynamicProbabilityVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

#endif //TRANSMISSION_NETWORKS_APP_MATRIX_H
