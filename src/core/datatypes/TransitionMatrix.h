//
// Created by Maxwell Murphy on 2/7/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSITIONMATRIX_H
#define TRANSMISSION_NETWORKS_APP_TRANSITIONMATRIX_H

#include "Matrix.h"

template <int MAX_COI>
using COITransitionMatrix = LogProbabilityTransitionMatrix<MAX_COI>;

template <int MAX_TRANSMISSIONS, int MAX_COI>
using TransmissionProbability = std::array<COITransitionMatrix<MAX_COI>, MAX_TRANSMISSIONS>;

#endif //TRANSMISSION_NETWORKS_APP_TRANSITIONMATRIX_H

