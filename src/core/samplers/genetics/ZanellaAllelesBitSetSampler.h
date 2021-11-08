//
// Created by Maxwell Murphy on 9/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/computation/PartialLikelihood.h"
#include "core/datatypes/Alleles.h"
#include "core/utils/numerics.h"


namespace transmission_nets::core::samplers::genetics {
    using Likelihood = core::computation::Likelihood;

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    class ZanellaAllelesBitSetSampler : public AbstractSampler {
    public:
        ZanellaAllelesBitSetSampler(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(const AllelesBitSetImpl& curr, std::vector<Likelihood>& neighborhoodLik) noexcept;

//        [[nodiscard]] double calcLogAcceptanceRatio(AllelesBitSetImpl& curr, AllelesBitSetImpl& prop) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:

        std::vector<Likelihood> calculateNeighborhoodLik() noexcept;

        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;

        boost::random::uniform_01<> uniform_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

    };

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::ZanellaAllelesBitSetSampler(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept :
            parameter_(parameter), target_(target), rng_(rng) {}

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    void ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::update() noexcept {
        const std::string stateId = "ZanellaAllelesBitSetSampler init";
        parameter_->saveState(stateId);

        auto curr = parameter_->value();
        auto currLik = target_->value();

        auto currNeighborhood = calculateNeighborhoodLik();
        auto currNeighborhoodSum = core::utils::logSumExp(currNeighborhood);

        if (currNeighborhoodSum <= -std::numeric_limits<Likelihood>::infinity()) {
            rejections_++;
            parameter_->restoreState(stateId);
        } else {
            auto proposal = sampleProposal(curr, currNeighborhood);
            parameter_->setValue(proposal);

            auto propLik = target_->value();

            if (propLik <= -std::numeric_limits<Likelihood>::infinity()) {
                rejections_++;
                parameter_->restoreState(stateId);
            } else {
                auto propNeighborhoodSum = core::utils::logSumExp(calculateNeighborhoodLik());

                const Likelihood logAcceptanceRatio = .5 * propLik + currNeighborhoodSum - .5 * currLik - propNeighborhoodSum;
                const Likelihood logProbAccept = log(uniform_dist_(*rng_));

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

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    double ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }



    template<typename T, typename Engine, typename AllelesBitSetImpl>
    AllelesBitSetImpl
    ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::sampleProposal(const AllelesBitSetImpl& curr, std::vector<Likelihood>& neighborhoodLik) noexcept {
        auto tmp = curr;
        auto normedNeighborhoodLik = core::utils::expNormalize(neighborhoodLik);
        auto unifSample = uniform_dist_(*rng_);

        int i = 0;
        double cumsum = 0;
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

        tmp.flip(i - 1);
        return tmp;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    std::vector<Likelihood> ZanellaAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::calculateNeighborhoodLik() noexcept {
        std::vector<Likelihood> neighborhood{};
        auto tmp = parameter_->value();
        neighborhood.reserve(tmp.totalAlleles());

        for (unsigned int i = 0; i < tmp.totalAlleles(); ++i) {
            tmp.flip(i);
            if (tmp.totalPositiveCount() > 0) {
                parameter_->saveState("flip1");
                parameter_->setValue(tmp);
                neighborhood.push_back(target_->value() * 0.5);
                if(std::isnan(neighborhood.back())) {
                    std::cerr << "Encountered NaN: " << target_->value() << std::endl;
                }
                parameter_->restoreState("flip1");
            } else {
                neighborhood.push_back(-std::numeric_limits<Likelihood>::infinity());
            }
            tmp.flip(i);
        }

        return neighborhood;
    }

}


#endif//TRANSMISSION_NETWORKS_APP_ZANELLAALLELESBITSETSAMPLER_H
