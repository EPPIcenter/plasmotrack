//
// Created by Maxwell Murphy on 6/11/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMADDEDGESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMADDEDGESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"
#include "core/parameters/TransmissionNetwork.h"
#include "core/samplers/AbstractSampler.h"

namespace transmission_nets::core::samplers::topology {

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    class RandomAddEdgeSampler : public AbstractSampler {
    public:
        RandomAddEdgeSampler(parameters::TransmissionNetwork<NodeValueImpl> &network, T &target, Engine *rng);

        void update() noexcept override;

//    ParentSet<NodeValueImpl> sampleProposal() noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

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

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeImpl>
    RandomAddEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeImpl>::RandomAddEdgeSampler(parameters::TransmissionNetwork<NodeImpl> &network, T &target, Engine *rng) :
            network_(network), target_(target), rng_(rng) {
        nodeIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, network_.nodes().size() - 1));
    }

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    void RandomAddEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::update() noexcept {
        const std::string stateId = "AddEdge1";
        Likelihood curLik = target_.value();

        auto parentNode = network_.nodes()[nodeIndexSamplingDist_(*rng_)];
        auto childNode = network_.nodes()[nodeIndexSamplingDist_(*rng_)];

        if (network_.parentSet(childNode)->value().size() < MAX_PARENT_SET_CARDINALITY and !network_.createsCycle(parentNode, childNode)) {
            network_.parentSet(childNode)->saveState(stateId);
            network_.addEdge(parentNode, childNode);

            const Likelihood acceptanceRatio = target_.value() - curLik;
            const Likelihood logProbAccept = log(uniformDist_(*rng_));
            const bool accept = logProbAccept <= acceptanceRatio;

            if (accept) {
                acceptances_++;
                network_.parentSet(childNode)->acceptState();
            } else {
                rejections_++;
                network_.parentSet(childNode)->restoreState(stateId);
                assert(curLik == target_.value());
            }


        } else {
            rejections_++;
        }
        total_updates_++;

    }

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomAddEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomAddEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::rejections() noexcept {
        return rejections_;
    }

    template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
    double RandomAddEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

}




#endif //TRANSMISSION_NETWORKS_APP_RANDOMADDEDGESAMPLER_H
