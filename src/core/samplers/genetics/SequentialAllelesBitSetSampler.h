//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SEQUENTIALALLELESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_SEQUENTIALALLELESAMPLER_H

#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/utils/generators/RandomSequence.h"

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include <memory>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    class SequentialAllelesBitSetSampler : public AbstractSampler {
    public:
        SequentialAllelesBitSetSampler(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(AllelesBitSetImpl curr, size_t idx) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> allele_index_sampling_dist_{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::SequentialAllelesBitSetSampler(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept : parameter_(parameter), target_(target), rng_(rng) {
        allele_index_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type(0, parameter_->value().totalAlleles() - 1));
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    void SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::update() noexcept {
        const auto rs = core::utils::generators::randomSequence(1, parameter_->value().totalAlleles(), rng_);
        for (const auto& i : rs) {
            const std::string stateId = "State1";
            Likelihood curLik         = target_->value();
            parameter_->saveState(stateId);
            const auto proposal = sampleProposal(parameter_->value(), i);

            assert(!target_->isDirty());
            parameter_->setValue(proposal);
            const auto acceptanceRatio = target_->value() - curLik;
            const auto logProbAccept   = log(uniform_dist_(*rng_));
            const bool accept          = logProbAccept <= acceptanceRatio;

            if (accept) {
                acceptances_++;
                parameter_->acceptState();
            } else {
                rejections_++;
                parameter_->restoreState(stateId);
            }

            assert(!target_->isDirty());

            total_updates_++;
        }
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    double SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    AllelesBitSetImpl
    SequentialAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::sampleProposal(const AllelesBitSetImpl curr, size_t idx) noexcept {
        auto tmp = curr;
//        tmp.flip(allele_index_sampling_dist_(*rng_));
//
//        while (tmp.totalPositiveCount() == 0) {
//            tmp.flip(allele_index_sampling_dist_(*rng_));
//        }
        tmp.flip(idx);

        if (tmp.totalPositiveCount() == 0) {
            tmp.flip(allele_index_sampling_dist_(*rng_));
        }

        return tmp;
    }
}// namespace transmission_nets::core::samplers::genetics


#endif//TRANSMISSION_NETWORKS_APP_SEQUENTIALALLELESAMPLER_H
