//
// A Run's chains: the Hot chain and its Tempered chains, each with its own State, Model and Sampler set.
//

#ifndef TRANSMISSION_NETWORKS_APP_ENSEMBLE_H
#define TRANSMISSION_NETWORKS_APP_ENSEMBLE_H

#include "core/samplers/meta/ReplicaExchange.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/SamplerSet.h"
#include "impl/model/Model/State.h"

namespace transmission_nets::impl::Model {
    using Ensemble = core::samplers::ReplicaExchange<State, Model, SamplerSet>;
}

#endif//TRANSMISSION_NETWORKS_APP_ENSEMBLE_H
