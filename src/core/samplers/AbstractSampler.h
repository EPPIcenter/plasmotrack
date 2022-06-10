//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H

#include "core/computation/PartialLikelihood.h"

namespace transmission_nets::core::samplers {
    using Likelihood = core::computation::Likelihood;
    class AbstractSampler {
    public:
        virtual ~AbstractSampler()     = default;
        virtual void update() noexcept = 0;
        virtual void adapt() noexcept {};
        virtual void adapt([[maybe_unused]] unsigned int idx) noexcept {};
    };

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
