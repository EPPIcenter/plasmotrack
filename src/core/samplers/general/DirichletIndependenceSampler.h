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

    namespace dirichlet {
        inline double logDensity(const std::vector<double>& alpha, const std::vector<double>& f) {
            double sumAlpha = 0, out = 0;
            for (size_t a = 0; a < alpha.size(); ++a) {
                sumAlpha += alpha[a];
                out += (alpha[a] - 1) * std::log(f[a]) - std::lgamma(alpha[a]);
            }
            return out + std::lgamma(sumAlpha);
        }

        // An equal mixture of Dirichlet(narrow) and Dirichlet(wide).
        inline double logMixtureDensity(const std::vector<double>& narrow, const std::vector<double>& wide, const std::vector<double>& f) {
            const double n = logDensity(narrow, f), w = logDensity(wide, f);
            const double m = std::max(n, w);
            return m + std::log(0.5 * std::exp(n - m) + 0.5 * std::exp(w - m));
        }

        // A draw from the mixture into `f`; false when an element falls below `lowerLimit`.
        template<typename Engine>
        bool drawMixture(const std::vector<double>& narrow, const std::vector<double>& wide, const double lowerLimit, Engine& rng, std::vector<double>& f) {
            const auto& alpha = boost::random::uniform_01<>()(rng) < 0.5 ? narrow : wide;
            f.resize(alpha.size());
            double total = 0;
            for (size_t a = 0; a < alpha.size(); ++a) {
                f[a] = boost::random::gamma_distribution<>(alpha[a])(rng);
                total += f[a];
            }
            bool inBounds = total > 0;
            for (auto& x : f) {
                x /= total;
                inBounds = inBounds and x >= lowerLimit;
            }
            return inBounds;
        }
    }// namespace dirichlet

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

    private:
        std::shared_ptr<parameters::Parameter<datatypes::Simplex>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        std::vector<double> alpha_;
        std::vector<double> wide_;
        double lowerLimit_;
        std::vector<double> draw_{};
    };

    template<typename T, typename Engine>
    void DirichletIndependenceSampler<T, Engine>::update() noexcept {
        auto proposal = this->propose(*target_, *rng_);
        if (not dirichlet::drawMixture(alpha_, wide_, lowerLimit_, *rng_, draw_)) {
            proposal.reject();
            return;
        }
        const auto current = parameter_->value().frequencies();
        proposal.change(*parameter_, datatypes::Simplex::exactly(draw_));
        proposal.metropolis(dirichlet::logMixtureDensity(alpha_, wide_, current) - dirichlet::logMixtureDensity(alpha_, wide_, draw_));
    }

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_DIRICHLETINDEPENDENCESAMPLER_H
