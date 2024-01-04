//
// Created by Maxwell Murphy on 4/24/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_DISCRETEDISTRIBUTION_H
#define TRANSMISSION_NETWORKS_APP_DISCRETEDISTRIBUTION_H

#include "core/parameters/Parameter.h"
#include "core/computation/PartialLikelihood.h"


namespace transmission_nets::core::distributions {

    using DiscreteDistribution = parameters::Parameter<std::vector<computation::Probability>>;

}// namespace transmission_nets::core::distributions

#endif//TRANSMISSION_NETWORKS_APP_DISCRETEDISTRIBUTION_H
