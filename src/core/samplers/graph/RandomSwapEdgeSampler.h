//
// Created by Maxwell Murphy on 6/11/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMSWAPEDGESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMSWAPEDGESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"

#include "core/samplers/AbstractSampler.h"

#include "core/containers/ParentSet.h"
#include "core/parameters/TransmissionNetwork.h"

namespace transmission_nets::core::samplers::graph {

    template<typename T, typename Engine, typename NodeValueImpl>
    class RandomSwapEdgeSampler : public AbstractSampler {
    public:
        RandomSwapEdgeSampler(parameters::TransmissionNetwork<NodeValueImpl> &network, T &target, Engine *rng);

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
        boost::random::uniform_int_distribution<> parentSetIndexSamplingDist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

    };

    template<typename T, typename Engine, typename NodeImpl>
    RandomSwapEdgeSampler<T, Engine, NodeImpl>::RandomSwapEdgeSampler(parameters::TransmissionNetwork<NodeImpl> &network, T &target, Engine *rng) :
            network_(network), target_(target), rng_(rng) {
        nodeIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, network_.nodes().size() - 1));
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    void RandomSwapEdgeSampler<T, Engine, NodeValueImpl>::update() noexcept {
        const std::string stateId = "SwapEdge1";
        double curLik = target_.value();

        const auto nodeA = network_.nodes()[nodeIndexSamplingDist_(*rng_)];
        const auto nodeB = network_.nodes()[nodeIndexSamplingDist_(*rng_)];
        if (nodeA != nodeB) {
            parameters::Parameter<containers::ParentSet<NodeValueImpl>>* nodeAParentSetParam = network_.parentSet(nodeA);
            parameters::Parameter<containers::ParentSet<NodeValueImpl>>* nodeBParentSetParam = network_.parentSet(nodeB);
            const auto& nodeAParentSet = nodeAParentSetParam->value();
            const auto& nodeBParentSet = nodeBParentSetParam->value();

            if (nodeAParentSet.size() > 0 and nodeBParentSet.size() > 0) {
                parentSetIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, nodeAParentSet.size() - 1));
                const unsigned int nodeAParentIdx = parentSetIndexSamplingDist_(*rng_);

                parentSetIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, nodeBParentSet.size() - 1));
                const unsigned int nodeBParentIdx = parentSetIndexSamplingDist_(*rng_);

                NodeValueImpl* nodeAParentNode = *(nodeAParentSet.begin() + nodeAParentIdx);
                NodeValueImpl* nodeBParentNode = *(nodeBParentSet.begin() + nodeBParentIdx);

                nodeAParentSetParam->saveState(stateId);
                assert(nodeAParentSetParam->isSaved());
                nodeBParentSetParam->saveState(stateId);
                assert(nodeBParentSetParam->isSaved());
                assert(!target_.isDirty());
                network_.removeEdge(nodeAParentNode, nodeA);
                network_.removeEdge(nodeBParentNode, nodeB);

                if (network_.createsCycle(nodeAParentNode, nodeB) or network_.createsCycle(nodeBParentNode, nodeA)) {
                    rejections_++;
                    nodeAParentSetParam->restoreState(stateId);
                    nodeBParentSetParam->restoreState(stateId);
                    assert(!target_.isDirty());
                    assert(curLik == target_.value());
                } else {
                    network_.addEdge(nodeAParentNode, nodeB);
                    network_.addEdge(nodeBParentNode, nodeA);
                    const double acceptanceRatio = target_.value() - curLik;
                    const double logProbAccept = log(uniformDist_(*rng_));
                    const bool accept = logProbAccept <= acceptanceRatio;
                    if (accept) {
                        acceptances_++;
                        nodeAParentSetParam->acceptState();
                        nodeBParentSetParam->acceptState();
                    } else {
                        rejections_++;
                        nodeAParentSetParam->restoreState(stateId);
                        nodeBParentSetParam->restoreState(stateId);
                        assert(!target_.isDirty());
                        assert(curLik == target_.value());
                    }
                }

            } else {
                rejections_++;
            }
        } else {
            rejections_++;
        }
        total_updates_++;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomSwapEdgeSampler<T, Engine, NodeValueImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomSwapEdgeSampler<T, Engine, NodeValueImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    double RandomSwapEdgeSampler<T, Engine, NodeValueImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

}


#endif //TRANSMISSION_NETWORKS_APP_RANDOMSWAPEDGESAMPLER_H
