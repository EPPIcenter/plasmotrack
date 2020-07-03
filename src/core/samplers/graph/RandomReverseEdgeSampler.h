//
// Created by Maxwell Murphy on 6/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMREVERSEEDGESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMREVERSEEDGESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"

#include "core/samplers/AbstractSampler.h"

#include "core/containers/ParentSet.h"
#include "core/containers/TransmissionNetwork.h"


template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
class RandomReverseEdgeSampler : public AbstractSampler {
public:
    RandomReverseEdgeSampler(TransmissionNetwork<NodeValueImpl> &network, T &target, Engine *rng);

    void update() noexcept override;

//    ParentSet<NodeValueImpl> sampleProposal() noexcept;

    [[nodiscard]] unsigned int acceptances() noexcept;

    [[nodiscard]] unsigned int rejections() noexcept;

    [[nodiscard]] double acceptanceRate() noexcept;

private:
    TransmissionNetwork<NodeValueImpl> &network_;
    T &target_;
    Engine *rng_;

    boost::random::uniform_01<> uniformDist_{};
    boost::random::uniform_int_distribution<> nodeIndexSamplingDist_;
    boost::random::uniform_int_distribution<> parentSetIndexSamplingDist_{};

    unsigned int acceptances_ = 0;
    unsigned int rejections_ = 0;
    unsigned int total_updates_ = 0;

};

template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeImpl>
RandomReverseEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeImpl>::RandomReverseEdgeSampler(TransmissionNetwork<NodeImpl> &network, T &target, Engine *rng) :
        network_(network), target_(target), rng_(rng) {
    nodeIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, network_.nodes().size() - 1));
}

template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
void RandomReverseEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::update() noexcept {
    double curLik = target_.value();

    auto childNode = network_.nodes()[nodeIndexSamplingDist_(*rng_)];
    auto childParentSetParam = network_.parentSet(childNode);
    const auto& childParentSet = childParentSetParam->value();
    if (childParentSet.size() != 0) {
        parentSetIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, childParentSet.size() - 1));

        const unsigned int idx = parentSetIndexSamplingDist_(*rng_);
        const auto& parentNode = *(childParentSet.begin() + idx);
        auto parentParentSetParam = network_.parentSet(parentNode);

        if (parentParentSetParam->value().size() < MAX_PARENT_SET_CARDINALITY) {
            childParentSetParam->saveState();
            parentParentSetParam->saveState();
            network_.removeEdge(parentNode, childNode);
            if (network_.createsCycle(childNode, parentNode)) {
                rejections_++;
                childParentSetParam->restoreState();
                parentParentSetParam->restoreState();
                assert(curLik == target_.value());
            } else {
                network_.addEdge(childNode, parentNode);
                const double acceptanceRatio = target_.value() - curLik;
                const double logProbAccept = log(uniformDist_(*rng_));
                const bool accept = logProbAccept <= acceptanceRatio;

                if (accept) {
                    acceptances_++;
                    childParentSetParam->acceptState();
                    parentParentSetParam->acceptState();
                } else {
                    rejections_++;
                    childParentSetParam->restoreState();
                    parentParentSetParam->restoreState();
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

template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
unsigned int RandomReverseEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::acceptances() noexcept {
    return acceptances_;
}

template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
unsigned int RandomReverseEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::rejections() noexcept {
    return rejections_;
}

template<int MAX_PARENT_SET_CARDINALITY, typename T, typename Engine, typename NodeValueImpl>
double RandomReverseEdgeSampler<MAX_PARENT_SET_CARDINALITY, T, Engine, NodeValueImpl>::acceptanceRate() noexcept {
    return double(acceptances_) / (acceptances_ + rejections_);
}

#endif //TRANSMISSION_NETWORKS_APP_RANDOMREVERSEEDGESAMPLER_H
