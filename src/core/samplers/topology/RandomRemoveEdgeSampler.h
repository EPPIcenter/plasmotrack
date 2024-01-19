//
// Created by Maxwell Murphy on 6/11/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMREMOVEEDGESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMREMOVEEDGESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"

#include "core/samplers/AbstractSampler.h"

#include "core/containers/ParentSet.h"
#include "core/parameters/TransmissionNetwork.h"

namespace transmission_nets::core::samplers::topology {

    template<typename T, typename Engine, typename NodeValueImpl>
    class RandomRemoveEdgeSampler : public AbstractSampler {
    public:
        RandomRemoveEdgeSampler(parameters::TransmissionNetwork<NodeValueImpl>& network, T& target, Engine* rng);

        void update() noexcept override;

        //    ParentSet<NodeValueImpl> sampleProposal() noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        parameters::TransmissionNetwork<NodeValueImpl>& network_;
        T& target_;
        Engine* rng_;

        boost::random::uniform_01<> uniformDist_{};
        boost::random::uniform_int_distribution<> nodeIndexSamplingDist_;
        boost::random::uniform_int_distribution<> parentSetIndexSamplingDist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine, typename NodeImpl>
    RandomRemoveEdgeSampler<T, Engine, NodeImpl>::RandomRemoveEdgeSampler(parameters::TransmissionNetwork<NodeImpl>& network, T& target, Engine* rng) : network_(network), target_(target), rng_(rng) {
        nodeIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, network_.nodes().size() - 1));
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    void RandomRemoveEdgeSampler<T, Engine, NodeValueImpl>::update() noexcept {
        const std::string stateId = "RemoveEdge1";
        Likelihood curLik         = target_.value();

        auto childNode        = network_.nodes()[nodeIndexSamplingDist_(*rng_)];
        auto param            = network_.parentSet(childNode);
        const auto& parentSet = param->value();
        if (parentSet.size() != 0) {
            parentSetIndexSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type(0, parentSet.size() - 1));

            const unsigned int idx = parentSetIndexSamplingDist_(*rng_);
            const auto& parentNode = *(parentSet.begin() + idx);

            param->saveState(stateId);
            network_.removeEdge(parentNode, childNode);

            const Likelihood acceptanceRatio = target_.value() - curLik;
            const Likelihood logProbAccept   = log(uniformDist_(*rng_));
            const bool accept                = logProbAccept <= acceptanceRatio;

            if (accept) {
                acceptances_++;
                param->acceptState();
            } else {
                rejections_++;
                param->restoreState(stateId);
                assert(curLik == target_.value());
            }

        } else {
            rejections_++;
        }
        total_updates_++;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomRemoveEdgeSampler<T, Engine, NodeValueImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    unsigned int RandomRemoveEdgeSampler<T, Engine, NodeValueImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename NodeValueImpl>
    double RandomRemoveEdgeSampler<T, Engine, NodeValueImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

}// namespace transmission_nets::core::samplers::topology


#endif//TRANSMISSION_NETWORKS_APP_RANDOMREMOVEEDGESAMPLER_H
