//
// The Run recipe's model choices: the leak prior, the error model, and ModelChoice, which carries them with the null
// model and the Transmission leak to the State. Each choice is parsed and given its canonical form here, and only here.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELCHOICE_H
#define TRANSMISSION_NETWORKS_APP_MODELCHOICE_H

#include <string>

namespace transmission_nets::impl::Model {
    // The prior on the Transmission leak ε when it is sampled. Parsed from "loguniform:LO,HI" (density 1/ε on [LO, HI];
    // "loguniform" alone is loguniform:1e-7,0.1) or "exponential:RATE" (truncated to (0, 0.5); "exponential" alone is
    // RATE 100).
    struct LeakPrior {
        static constexpr double EXPONENTIAL_UPPER = 0.5;

        enum class Kind { LogUniform, Exponential } kind = Kind::LogUniform;
        // LogUniform's range.
        double lo = 1e-7;
        double hi = 0.1;
        // Exponential's rate.
        double rate = 100;

        // The range the leak's walk keeps it strictly inside.
        [[nodiscard]] double lower() const { return kind == Kind::LogUniform ? lo : 0.0; }
        [[nodiscard]] double upper() const { return kind == Kind::LogUniform ? hi : EXPONENTIAL_UPPER; }
        [[nodiscard]] double logDensity(double leak) const;
        // Whether the walk can start at `leak`: strictly between lower() and upper().
        [[nodiscard]] bool holds(double leak) const { return leak > lower() and leak < upper(); }
        // `leak` if the prior holds it; otherwise the nearest value just inside the range.
        [[nodiscard]] double startAt(double leak) const;
    };
    // Throws std::invalid_argument for a description it can't read, a number with anything after it included.
    LeakPrior parseLeakPrior(const std::string& description);
    // The description in one form for each prior ("loguniform" and "loguniform:1e-07,0.1" give the same), or the
    // description unchanged if it can't be read.
    std::string canonicalLeakPrior(const std::string& description);

    // How the observation error rates are shared between Infections:
    //   infection    each Infection its own two rates, prior Gamma(10, 0.001) (the original model);
    //   shared       one false-positive and one false-negative rate for every Infection, priors Gamma(2, 0.01) and
    //                Gamma(2, 0.02) (the WIDE_ constants in config.h);
    //   hierarchical each Infection its own rates, drawn from Gamma(shape, mean / shape) with the mean and shape of
    //                each rate sampled (Run scalars efp_mean, efp_shape, efn_mean, efn_shape).
    enum class ErrorModel { Infection, Shared, Hierarchical };
    // Throws std::invalid_argument for a name it doesn't know; "" is `infection`.
    ErrorModel parseErrorModel(const std::string& name);
    // The name in one form for each model ("" for `infection`, given or not), or the name unchanged if it isn't one.
    std::string canonicalErrorModel(const std::string& name);

    // The Run recipe's choices that shape each chain's State and Sampler set. Not Quantities: never sampled or
    // checkpointed.
    struct ModelChoice {
        bool nullModel = false;
        // The Transmission leak ε (see doc/marginal-latent-parent.md, section 6); with a leak prior, its start.
        double transmissionLeak = 0;
        // Empty: the leak is fixed at transmissionLeak. Otherwise a LeakPrior description: the leak is a Run scalar
        // sampled under it, starting at LeakPrior::startAt(transmissionLeak).
        std::string leakPrior{};
        // An ErrorModel name; empty is `infection`.
        std::string errorModel{};
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_MODELCHOICE_H
