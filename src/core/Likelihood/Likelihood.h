//
// Created by Maxwell Murphy on 12/9/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_LIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_LIKELIHOOD_H

#include "vector"
#include "../Node.h"
#include "PartialLikelihood.h"

class Likelihood : public Node<float> {
private:
    boost::container::flat_set<PartialLikelihood *> likelihoods_;
    boost::container::flat_set<PartialLikelihood *> dirty_dependencies_{}; // Dirty elements this node depends on

public:
    Likelihood();

    void pushLikelihood(PartialLikelihood *ptr);

    float value() override;

    float peek() override;

    void reinitialize();

};


#endif //TRANSMISSION_NETWORKS_APP_LIKELIHOOD_H
