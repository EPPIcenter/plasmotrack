//
// Created by Maxwell Murphy on 10/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEXSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEXSAMPLER_H

#include <boost/random.hpp>

#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

namespace transmission_nets::core::samplers {

    template<typename T, typename Engine = boost::random::mt19937>
    class SimplexSampler : public AbstractSampler {

    public:
        SimplexSampler(parameters::Parameter<datatypes::Simplex>& parameter, T& target, Engine* rng);
        SimplexSampler(parameters::Parameter<datatypes::Simplex>& parameter, T& target, Engine* rng, double variance);

        void setAdaptationRate(double adaptationRate);

        void setTargetAcceptanceRate(double targetAcceptanceRate);

        void setMinVariance(double minVariance);

        void setMaxVariance(double maxVariance);

        void update() noexcept override;

        void adapt() noexcept override;

        void adapt(unsigned int idx) noexcept override;

        double acceptanceRate(int idx) noexcept;

    private:
        parameters::Parameter<datatypes::Simplex>& parameter_;
        T& target_;
        Engine* rng_;
        boost::random::normal_distribution<> normal_dist_{0, 1};
        boost::random::uniform_01<> uniform_dist_{};

        double variance_    = 1;
        double acceptances_ = 0;
        double rejections_  = 0;

        double minVariance_ = 1e-12;
        double maxVariance_ = 1e6;

        double adaptation_rate_        = 1;
        double target_acceptance_rate_ = .23;

        int total_updates_ = 0;
    };

    template<typename T, typename Engine>
    SimplexSampler<T, Engine>::SimplexSampler(parameters::Parameter<datatypes::Simplex>& parameter, T& target, Engine* rng) : parameter_(parameter), target_(target), rng_(rng) {}

    template<typename T, typename Engine>
    SimplexSampler<T, Engine>::SimplexSampler(parameters::Parameter<datatypes::Simplex>& parameter, T& target, Engine* rng, double variance) : parameter_(parameter), target_(target), rng_(rng), variance_(variance) {}

    template<typename T, typename Engine>
    void SimplexSampler<T, Engine>::setAdaptationRate(double adaptationRate) {
        adaptation_rate_ = adaptationRate;
    }


}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SIMPLEXSAMPLER_H
