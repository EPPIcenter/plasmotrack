//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/datatypes/Alleles.h"


namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    class RandomAllelesBitSetSampler : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler(std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        void update() noexcept override;

        AllelesBitSetImpl sampleProposal(AllelesBitSetImpl curr) noexcept;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> allele_index_sampling_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

    };

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::RandomAllelesBitSetSampler(
            std::shared_ptr<parameters::Parameter<AllelesBitSetImpl>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept :
            parameter_(parameter), target_(target), rng_(rng) {
        allele_index_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type(0, parameter_->value().totalAlleles() - 1)
        );
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    void RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::update() noexcept {
        const std::string stateId = "State1";
        Likelihood curLik = target_->value();
        parameter_->saveState(stateId);
        const auto proposal = sampleProposal(parameter_->value());

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const auto acceptanceRatio = target_->value() - curLik;
        const auto logProbAccept = log(uniform_dist_(*rng_));
        const bool accept = logProbAccept <= acceptanceRatio;

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

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    unsigned int RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    double RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::acceptanceRate() noexcept {
        return double(acceptances_) / (acceptances_ + rejections_);
    }

    template<typename T, typename Engine, typename AllelesBitSetImpl>
    AllelesBitSetImpl
    RandomAllelesBitSetSampler<T, Engine, AllelesBitSetImpl>::sampleProposal(const AllelesBitSetImpl curr) noexcept {
        auto tmp = curr;
        tmp.flip(allele_index_sampling_dist_(*rng_));

        while(tmp.totalPositiveCount() == 0) {
            tmp.flip(allele_index_sampling_dist_(*rng_));
        }

        return tmp;
    }
}



#endif //TRANSMISSION_NETWORKS_APP_ALLELESAMPLER_H
