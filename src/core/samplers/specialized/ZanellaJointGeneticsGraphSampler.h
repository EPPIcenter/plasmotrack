//
// Created by Maxwell Murphy on 9/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZANELLAJOINTGENETICSGRAPHSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ZANELLAJOINTGENETICSGRAPHSAMPLER_H

#include <boost/random.hpp>

#include "core/parameters/TransmissionNetwork.h"
#include "core/samplers/AbstractSampler.h"

namespace transmission_nets::core::samplers {

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    class ZanellaJointGeneticsGraphSampler : public AbstractSampler {
    public:
    private:

        parameters::TransmissionNetwork<NodeValueImpl> &network_;
        T &target_;
        Engine *rng_;

        boost::random::uniform_01<> uniformDist_{};
        boost::random::uniform_int_distribution<> nodeIndexSamplingDist_;

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

    };

}

#endif//TRANSMISSION_NETWORKS_APP_ZANELLAJOINTGENETICSGRAPHSAMPLER_H
