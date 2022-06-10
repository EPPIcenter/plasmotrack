//
// Created by Maxwell Murphy on 10/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZANELLAORDERSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ZANELLAORDERSAMPLER_H

#include <boost/random.hpp>

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Ordering.h"
#include "core/samplers/AbstractSampler.h"
#include "core/utils/numerics.h"

namespace transmission_nets::core::samplers {
    using Likelihood = core::computation::Likelihood;
    /*
     * Scheduler using a locally informed transition kernel. Randomly selects an index then samples from swaps no more than k indices away.
     */

    template<typename T, typename OrderingElement, typename Engine = boost::random::mt19937>
    class ZanellaOrderSampler : public AbstractSampler {
    public:
        ZanellaOrderSampler(parameters::Ordering<OrderingElement>& parameter, T& target, Engine* rng) noexcept;

        ZanellaOrderSampler(parameters::Ordering<OrderingElement>& parameter, T& target, Engine* rng, int neighborhoodSize) noexcept;

        void update() noexcept override;

        std::pair<std::vector<int>, std::vector<Likelihood>> calculateNeighborhoodLik(int elementIdx);

        int sampleProposal(const std::vector<int>& elements, const std::vector<Likelihood>& neighborhoodLik) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

        void setNeighborhoodSize(int neighborhoodSize) noexcept;

    private:
        parameters::Ordering<OrderingElement>& parameter_;
        T& target_;
        Engine* rng_;
        int numElements_;
        int neighborhoodSize_ = 4;
        boost::random::uniform_int_distribution<> idxSamplingDist_{};
        boost::random::uniform_01<> uniformDist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };


    template<typename T, typename OrderingElement, typename Engine>
    void ZanellaOrderSampler<T, OrderingElement, Engine>::update() noexcept {
        const std::string stateId = "ZanellaOrderSampler";
        parameter_.saveState(stateId);

        const auto curr          = parameter_.value();
        const Likelihood currLik = target_.value();
        int elementIdx           = idxSamplingDist_(*rng_);

        auto [elements, currNeighborhood]       = calculateNeighborhoodLik(elementIdx);
        const Likelihood currNeighborhoodLogSum = core::utils::logSumExp(currNeighborhood);

        if (currNeighborhoodLogSum <= -std::numeric_limits<Likelihood>::infinity()) {
            rejections_++;
            parameter_.restoreState(stateId);
        } else {
            auto proposal = sampleProposal(elements, currNeighborhood);
            parameter_.swap(elementIdx, proposal);

            const Likelihood propLik = target_.value();

            if (propLik <= -std::numeric_limits<Likelihood>::infinity()) {
                rejections_++;
                parameter_.restoreState(stateId);
            } else {
                auto [proposalElements, proposalNeighborhood] = calculateNeighborhoodLik(proposal);
                const Likelihood propNeighborhoodLogSum       = core::utils::logSumExp(proposalNeighborhood);

                const Likelihood logAcceptanceRatio = .5 * propLik + currNeighborhoodLogSum - .5 * currLik - propNeighborhoodLogSum;
                const Likelihood logProbAccept      = log(uniformDist_(*rng_));

                const bool accept = (logProbAccept <= logAcceptanceRatio);

                if (accept) {
                    acceptances_++;
                    parameter_.acceptState();
                } else {
                    rejections_++;
                    parameter_.restoreState(stateId);
                    assert(currLik == target_.value());
                }
            }
        }


        assert(!target_.isDirty());

        total_updates_++;
    }


    template<typename T, typename OrderingElement, typename Engine>
    ZanellaOrderSampler<T, OrderingElement, Engine>::ZanellaOrderSampler(parameters::Ordering<OrderingElement>& parameter, T& target, Engine* rng) noexcept
        : parameter_(parameter), target_(target), rng_(rng) {
        numElements_ = parameter.value().size();
        idxSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type{0, static_cast<int>(numElements_ - 1)});
    }


    template<typename T, typename OrderingElement, typename Engine>
    ZanellaOrderSampler<T, OrderingElement, Engine>::ZanellaOrderSampler(parameters::Ordering<OrderingElement>& parameter, T& target, Engine* rng, int neighborhoodSize) noexcept
        : parameter_(parameter), target_(target), rng_(rng), neighborhoodSize_(neighborhoodSize) {
        numElements_ = parameter.value().size();
        idxSamplingDist_.param(boost::random::uniform_int_distribution<>::param_type{0, static_cast<int>(numElements_ - 1)});
    }


    template<typename T, typename OrderingElement, typename Engine>
    int ZanellaOrderSampler<T, OrderingElement, Engine>::sampleProposal(const std::vector<int>& elements, const std::vector<Likelihood>& Lliks) noexcept {
        auto normedLlik = core::utils::expNormalize(Lliks);
        auto unifSample = uniformDist_(*rng_);

        int i             = 0;
        Likelihood cumsum = 0;
        while (unifSample >= cumsum) {
            if (i == int(normedLlik.size())) {
                std::cerr << "CumSum: " << cumsum << " " << unifSample << std::endl;
                std::cerr << "[";
                for (const auto llik : normedLlik) {
                    std::cerr << llik << ", ";
                }
                std::cerr << "]" << std::endl;

                std::cerr << "[";
                for (const auto llik : Lliks) {
                    std::cerr << llik << ", ";
                }
                std::cerr << "]" << std::endl;
            }
            assert(i < int(normedLlik.size()));
            cumsum += normedLlik[i];
            ++i;
        }

        return elements[i - 1];
    }


    template<typename T, typename OrderingElement, typename Engine>
    unsigned int ZanellaOrderSampler<T, OrderingElement, Engine>::acceptances() noexcept {
        return acceptances_;
    }


    template<typename T, typename OrderingElement, typename Engine>
    unsigned int ZanellaOrderSampler<T, OrderingElement, Engine>::rejections() noexcept {
        return rejections_;
    }


    template<typename T, typename OrderingElement, typename Engine>
    double ZanellaOrderSampler<T, OrderingElement, Engine>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }


    template<typename T, typename OrderingElement, typename Engine>
    std::pair<std::vector<int>, std::vector<Likelihood>> ZanellaOrderSampler<T, OrderingElement, Engine>::calculateNeighborhoodLik(int elementIdx) {

        std::vector<int> elements{};
        std::vector<Likelihood> neighborhood{};
        neighborhood.reserve(neighborhoodSize_ * 2);
        elements.reserve(neighborhoodSize_ * 2);

        int i     = std::max(0, elementIdx - neighborhoodSize_);
        int max_i = std::min(numElements_, elementIdx + neighborhoodSize_ + 1);

        for (; i < max_i; ++i) {
            if (i != elementIdx) {
                parameter_.saveState("swap");
                parameter_.swap(elementIdx, i);
                neighborhood.push_back(target_.value() * 0.5);
                if (std::isnan(neighborhood.back())) {
                    std::cerr << "Encountered NaN: " << target_.value() << std::endl;
                }
                elements.push_back(i);
                parameter_.restoreState("swap");
            }
        }

        return std::make_pair(elements, neighborhood);
    }


    template<typename T, typename OrderingElement, typename Engine>
    void ZanellaOrderSampler<T, OrderingElement, Engine>::setNeighborhoodSize(int neighborhoodSize) noexcept {
        neighborhoodSize_ = neighborhoodSize;
    }

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_ZANELLAORDERSAMPLER_H
