//
// Created by Maxwell Murphy on 9/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/computation/PartialLikelihood.h"
#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"


namespace transmission_nets::core::samplers::genetics {
    using Likelihood = core::computation::Likelihood;

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize = 2>
    class ZanellaAllelesBitSetSampler : public AbstractSampler {
    public:
        ZanellaAllelesBitSetSampler(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(const AllelesBitSetImpl& curr, const std::vector<Likelihood>& neighborhoodLik) noexcept;
        AllelesBitSetImpl sampleProposal(const AllelesBitSetImpl& curr, const std::vector<Likelihood>& neighborhoodLik, const std::vector<std::vector<unsigned int>>& indices) noexcept;

        //        [[nodiscard]] double calcLogAcceptanceRatio(AllelesBitSetImpl& curr, AllelesBitSetImpl& prop) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        std::vector<Likelihood> calculateNeighborhoodLik() noexcept;
        std::tuple<std::vector<std::vector<unsigned int>>, std::vector<Likelihood>> calculateNeighborhoodLik(unsigned int totalFlips) noexcept;

        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;

        boost::random::uniform_01<> uniform_dist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::ZanellaAllelesBitSetSampler(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept : parameter_(parameter), target_(target), rng_(rng) {}

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    void ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::update() noexcept {
        SAMPLER_STATE_ID stateId = SAMPLER_STATE_ID::ZanellaAllelesBitSetID;
        parameter_->saveState(stateId);

        const auto curr    = parameter_->value();
        const auto currLik = target_->value();

        auto currNeighborhood    = calculateNeighborhoodLik();
//        auto [indices, currNeighborhood] = calculateNeighborhoodLik(std::min(curr.totalAlleles() - 1, NeighborhoodSize));
        auto currNeighborhoodSum = core::utils::logSumExp(currNeighborhood);

        if (currNeighborhoodSum <= -std::numeric_limits<Likelihood>::infinity()) {
            rejections_++;
            parameter_->restoreState(stateId);
        } else {
            auto proposal = sampleProposal(curr, currNeighborhood);
//            auto proposal = sampleProposal(curr, currNeighborhood, indices);
            parameter_->setValue(proposal);

            auto propLik = target_->value();

            if (propLik <= -std::numeric_limits<Likelihood>::infinity()) {
                rejections_++;
                parameter_->restoreState(stateId);
            } else {
                auto propNeighborhoodSum = core::utils::logSumExp(calculateNeighborhoodLik());

                const Likelihood logAcceptanceRatio = .5 * propLik + currNeighborhoodSum - .5 * currLik - propNeighborhoodSum;
                const Likelihood logProbAccept      = log(uniform_dist_(*rng_));

                const bool accept = (logProbAccept <= logAcceptanceRatio);

                if (accept) {
                    acceptances_++;
                    parameter_->acceptState();
                } else {
                    rejections_++;
                    parameter_->restoreState(stateId);
                }
            }
        }


        assert(!target_->isDirty());
        total_updates_++;
        //        std::cout << "Acceptance Rate: " << acceptanceRate() << std::endl;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    unsigned int ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    unsigned int ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    double ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }


    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    AllelesBitSetImpl
    ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::sampleProposal(const AllelesBitSetImpl& curr, const std::vector<Likelihood>& neighborhoodLik) noexcept {
        auto tmp                   = curr;
        auto normedNeighborhoodLik = core::utils::expNormalize(neighborhoodLik);
        auto unifSample            = uniform_dist_(*rng_);

        int i         = 0;
        double cumsum = 0;
        while (unifSample >= cumsum) {

#ifdef DNDEBUG
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
#endif
            assert(i < int(normedNeighborhoodLik.size()));
            cumsum += normedNeighborhoodLik[i];
            ++i;
        }

        tmp.flip(i - 1);
        return tmp;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    AllelesBitSetImpl
    ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::sampleProposal(const AllelesBitSetImpl& curr, const std::vector<Likelihood>& neighborhoodLik, const std::vector<std::vector<unsigned int>>& indices) noexcept {
        auto tmp                   = curr;
        auto normedNeighborhoodLik = core::utils::expNormalize(neighborhoodLik);
        auto unifSample            = uniform_dist_(*rng_);

        int i         = 0;
        double cumsum = 0;
        while (unifSample >= cumsum) {

#ifdef DNDEBUG
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
#endif
            assert(i < int(normedNeighborhoodLik.size()));
            cumsum += normedNeighborhoodLik[i];
            ++i;
        }

        tmp.flip(indices[i - 1]);
        return tmp;
    }


    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    std::vector<Likelihood> ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::calculateNeighborhoodLik() noexcept {
        std::vector<Likelihood> neighborhood{};
        auto tmp = parameter_->value();
        neighborhood.reserve(tmp.totalAlleles());

        for (unsigned int i = 0; i < tmp.totalAlleles(); ++i) {
            tmp.flip(i);
            if (tmp.totalPositiveCount() > 0) {
                parameter_->saveState(SAMPLER_STATE_ID::ZanellaAllelesBitSetFlipID);
                parameter_->setValue(tmp);
                neighborhood.push_back(target_->value() * 0.5);
#ifdef DNDEBUG
                if (std::isnan(neighborhood.back())) {
                    std::cerr << "Encountered NaN: " << target_->value() << std::endl;
                }
#endif
                parameter_->restoreState(SAMPLER_STATE_ID::ZanellaAllelesBitSetFlipID);
            } else {
                neighborhood.push_back(-std::numeric_limits<Likelihood>::infinity());
            }
            tmp.flip(i);
        }

        return neighborhood;
    }


    template<typename T, typename Engine, typename AllelesBitSetImpl, unsigned int NeighborhoodSize>
    std::tuple<std::vector<std::vector<unsigned int>>, std::vector<Likelihood>> ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl, NeighborhoodSize>::calculateNeighborhoodLik(unsigned int totalFlips) noexcept {
        std::vector<Likelihood> neighborhood{};
        std::vector<std::vector<unsigned int>> combos{};
        auto tmp = parameter_->value();
        core::utils::generators::CombinationIndicesGenerator comboGenerator;

        for (unsigned int i = 1; i <= totalFlips; ++i) {
            comboGenerator.reset(tmp.totalAlleles(), i);
            while (!comboGenerator.completed) {
                auto combo = comboGenerator.curr;
                tmp.flip(combo);
                if (tmp.totalPositiveCount() > 0) {
                    parameter_->saveState("flip1");
                    parameter_->setValue(tmp);
                    neighborhood.push_back(target_->value() * 0.5);
                    combos.push_back(combo);
                    parameter_->restoreState("flip1");
                }
                tmp.flip(combo);
                comboGenerator.next();
            }
        }
        return std::make_tuple(combos, neighborhood);
    }

}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H
