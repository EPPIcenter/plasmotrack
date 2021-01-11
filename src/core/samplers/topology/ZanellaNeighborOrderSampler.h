//
// Created by Maxwell Murphy on 9/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZANELLANEIGHBORORDERSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ZANELLANEIGHBORORDERSAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

#include "core/utils/numerics.h"
#include "core/samplers/AbstractSampler.h"
#include "core/parameters/Ordering.h"
#include "core/io/serialize.h"

namespace transmission_nets::core::samplers {

    /*
     * Sampler using a locally informed transition kernel. Samples from the set of swapped neighboring elements in a given ordering.
     */
    template<typename T,  typename OrderingElement, typename Engine=boost::random::mt19937>
    class ZanellaNeighborOrderSampler : public AbstractSampler {
    public:
        ZanellaNeighborOrderSampler(parameters::Ordering<OrderingElement> &parameter, T &target, Engine *rng) noexcept;

        void update() noexcept override;

        std::vector<Likelihood> calculateNeighborhoodLik();

        int sampleProposal(const std::vector<Likelihood>& neighborhoodLik) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        parameters::Ordering<OrderingElement> &parameter_;
        T &target_;
        Engine *rng_;
        unsigned int num_elements_;
        boost::random::uniform_01<> uniform_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;
    };


    template<typename T, typename OrderingElement, typename Engine>
    void ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::update() noexcept {
        const std::string stateId = "ZanellaOrderSampler";
        parameter_.saveState(stateId);

        const auto curr = parameter_.value();
        const Likelihood currLik = target_.value();

        const auto currNeighborhood = calculateNeighborhoodLik();
        const auto currNeighborhoodLogSum = core::utils::logSumExp(currNeighborhood);

        if (currNeighborhoodLogSum <= -std::numeric_limits<Likelihood>::infinity()) {
            rejections_++;
            parameter_.restoreState(stateId);
        } else {
            auto proposal = sampleProposal(currNeighborhood);
            parameter_.swap(proposal, proposal + 1);

            const Likelihood propLik = target_.value();

            if (propLik <= -std::numeric_limits<Likelihood>::infinity()) {
                rejections_++;
                parameter_.restoreState(stateId);
            } else {
                const Likelihood propNeighborhoodLogSum = core::utils::logSumExp(calculateNeighborhoodLik());

                const Likelihood logAcceptanceRatio = .5 * propLik + currNeighborhoodLogSum - .5 * currLik - propNeighborhoodLogSum;
                const Likelihood logProbAccept = log(uniform_dist_(*rng_));

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
    ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::ZanellaNeighborOrderSampler(parameters::Ordering<OrderingElement> &parameter, T &target, Engine *rng) noexcept
            :parameter_(parameter), target_(target), rng_(rng) {
        num_elements_ = parameter.value().size();
    }


    template<typename T, typename OrderingElement, typename Engine>
    int ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::sampleProposal(const std::vector<Likelihood>& neighborhoodLik) noexcept {
        const auto normedNeighborhoodLik = core::utils::expNormalize(neighborhoodLik);
        auto unifSample = uniform_dist_(*rng_);

        int i = 0;
        Likelihood cumsum = 0;
        while (unifSample >= cumsum) {
            if (i == int(normedNeighborhoodLik.size())) {
                std::cerr << "CumSum: " << cumsum << " " << unifSample << std::endl;
                std::cerr << "[";
                for (const auto llik : normedNeighborhoodLik) {
                    std::cerr << llik << ", ";
                }
                std::cerr << "]" << std::endl;

                std::cerr << "[";
                for (const auto llik : neighborhoodLik) {
                    std::cerr << llik << ", ";
                }
                std::cerr << "]" << std::endl;
            }
            assert(i < int(normedNeighborhoodLik.size()));
            cumsum += normedNeighborhoodLik[i];
            ++i;
        }

        return i - 1;
    }

    template<typename T, typename OrderingElement, typename Engine>
    unsigned int ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename OrderingElement, typename Engine>
    unsigned int ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename OrderingElement, typename Engine>
    double ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

    template<typename T, typename OrderingElement, typename Engine>
    std::vector<Likelihood> ZanellaNeighborOrderSampler<T, OrderingElement, Engine>::calculateNeighborhoodLik() {

        std::vector<Likelihood> neighborhood{};
        neighborhood.reserve(num_elements_ - 1);
        for (size_t i = 0; i < parameter_.value().size() - 1; ++i) {
            parameter_.saveState("swap");
            parameter_.swap(i, i + 1);
            neighborhood.push_back(target_.value() * 0.5);
            if(std::isnan(neighborhood.back())) {
                std::cerr << "Encountered NaN: " << target_.value() << std::endl;
            }
            parameter_.restoreState("swap");
        }

        return neighborhood;
    }

}

#endif//TRANSMISSION_NETWORKS_APP_ZANELLANEIGHBORORDERSAMPLER_H
