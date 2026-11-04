//
// Proposes a whole simplex at once, independently of its current value, from a fixed Dirichlet mixture.
//

#ifndef TRANSMISSION_NETWORKS_APP_DIRICHLETINDEPENDENCESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_DIRICHLETINDEPENDENCESAMPLER_H

#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

#include <boost/random.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers {

    /*
     * An independence Metropolis-Hastings move on a simplex. The proposal is an equal mixture of Dirichlet(alpha) and
     * Dirichlet(widen * alpha), the same mean with a heavier spread, so that a value in alpha's tails is still proposed
     * often enough to be left. With alpha fitted to the parameter's posterior, an accepted move is close to an
     * independent draw, where a random walk moves one coordinate a little. A proposal with any element below the lower
     * limit is rejected, as the target has no mass there.
     */
    template<typename T, typename Engine = boost::random::mt19937>
    class DirichletIndependenceSampler : public AbstractSampler {
    public:
        DirichletIndependenceSampler(std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter,
                                     std::shared_ptr<T> target,
                                     std::shared_ptr<Engine> rng,
                                     std::vector<double> alpha,
                                     double widen,
                                     double lowerLimit)
            : parameter_(std::move(parameter)), target_(std::move(target)), rng_(std::move(rng)), alpha_(std::move(alpha)), wide_(alpha_), lowerLimit_(lowerLimit) {
            if (alpha_.size() != parameter_->value().totalElements()) {
                throw std::invalid_argument("Dirichlet concentrations do not match the simplex");
            }
            for (auto& a : wide_) {
                a *= widen;
            }
        }

        void update() noexcept override;

        [[nodiscard]] double logProposalDensity(const std::vector<double>& f) const;

    private:
        static double logDirichlet(const std::vector<double>& alpha, const std::vector<double>& f);

        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        std::vector<double> alpha_;
        std::vector<double> wide_;
        double lowerLimit_;
        boost::random::uniform_01<> uniform_{};
    };

    template<typename T, typename Engine>
    double DirichletIndependenceSampler<T, Engine>::logDirichlet(const std::vector<double>& alpha, const std::vector<double>& f) {
        double sumAlpha = 0, out = 0;
        for (size_t a = 0; a < alpha.size(); ++a) {
            sumAlpha += alpha[a];
            out += (alpha[a] - 1) * std::log(f[a]) - std::lgamma(alpha[a]);
        }
        return out + std::lgamma(sumAlpha);
    }

    template<typename T, typename Engine>
    double DirichletIndependenceSampler<T, Engine>::logProposalDensity(const std::vector<double>& f) const {
        const double narrow = logDirichlet(alpha_, f), wide = logDirichlet(wide_, f);
        const double m      = std::max(narrow, wide);
        return m + std::log(0.5 * std::exp(narrow - m) + 0.5 * std::exp(wide - m));
    }

    template<typename T, typename Engine>
    void DirichletIndependenceSampler<T, Engine>::update() noexcept {
        const auto& alpha = uniform_(*rng_) < 0.5 ? alpha_ : wide_;
        std::vector<double> f(alpha.size());
        double total = 0;
        for (size_t a = 0; a < alpha.size(); ++a) {
            f[a] = boost::random::gamma_distribution<>(alpha[a])(*rng_);
            total += f[a];
        }
        bool outOfBounds = not(total > 0);
        for (auto& x : f) {
            x /= total;
            outOfBounds = outOfBounds or not(x >= lowerLimit_);
        }
        auto proposal = this->propose(*target_, *rng_);
        if (outOfBounds) {
            proposal.reject();
            return;
        }
        const auto current = parameter_->value().frequencies();
        proposal.change(*parameter_, datatypes::Simplex::exactly(f));
        proposal.metropolis(logProposalDensity(current) - logProposalDensity(f));
    }

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_DIRICHLETINDEPENDENCESAMPLER_H
