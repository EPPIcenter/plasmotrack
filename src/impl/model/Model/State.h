//
// The State: every Quantity a chain samples, built from a Study. It applies the model's rules to the data (starting
// genotypes capped at the COI limit, the frequency floor, the null model) and refuses an input that breaks one.
// quantities(visitor) enumerates the Quantities once; the Checkpoint and the Run recorder's Traces follow that
// enumeration, and the Sampler set attaches its samplers by the same names.
//

#ifndef TRANSMISSION_NETWORKS_APP_STATE_H
#define TRANSMISSION_NETWORKS_APP_STATE_H

#include "config.h"
#include "Refusal.h"
#include "Study.h"
#include "core/containers/AllowedRelationships.h"
#include "core/distributions/DiscreteDistribution.h"

#include <boost/random.hpp>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>


namespace transmission_nets::impl::Model {
    // The prior on the transmission leak ε when it is sampled. Parsed from "loguniform:LO,HI" (density 1/ε on [LO, HI];
    // "loguniform" alone is loguniform:1e-7,0.1) or "exponential:RATE" (on (0, 0.5); "exponential" alone is RATE 100).
    struct LeakPrior {
        enum class Kind { LogUniform, Exponential } kind = Kind::LogUniform;
        double a = 1e-7;// LO, or RATE
        double b = 0.1; // HI
        // The range the sampler walks in.
        [[nodiscard]] double lower() const { return kind == Kind::LogUniform ? a : 0.0; }
        [[nodiscard]] double upper() const { return kind == Kind::LogUniform ? b : 0.5; }
        [[nodiscard]] double logDensity(double leak) const;
    };
    // Throws std::invalid_argument for a description it can't read.
    LeakPrior parseLeakPrior(const std::string& description);
    // The description in one form for each prior ("loguniform" and "loguniform:1e-07,0.1" give the same), or the
    // description unchanged if it can't be read.
    std::string canonicalLeakPrior(const std::string& description);

    // How the observation error rates are shared between Infections:
    //   infection    each Infection its own two rates, prior Gamma(10, 0.001) (the original model);
    //   shared       one false-positive and one false-negative rate for every Infection, priors Gamma(2, 0.01) and
    //                Gamma(2, 0.02);
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
        // The Transmission process's leak ε (see doc/marginal-latent-parent.md, section 6).
        double transmissionLeak = 0;
        // Empty: the leak is fixed at transmissionLeak. Otherwise a LeakPrior description: the leak is a Run scalar
        // sampled under it, starting at transmissionLeak (moved inside the prior's range).
        std::string leakPrior{};
        // An ErrorModel name; empty is `infection`.
        std::string errorModel{};
    };

    struct State {
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        using p_DiscreteDist = std::shared_ptr<core::distributions::DiscreteDistribution>;

        // Builds the starting State of `study`: each unobserved Locus's Latent genotype is one random Allele drawn
        // from `rng`. Throws Refused for a Study whose Infections name a Locus or a parent it doesn't have, or whose
        // genotypes or frequencies don't fit their Locus.
        State(const Study& study, const std::shared_ptr<boost::random::mt19937>& rng, const ModelChoice& choice);
        // The same, for these choices and the rest at their defaults.
        explicit State(const Study& study, const std::shared_ptr<boost::random::mt19937>& rng, bool null_model = false, double transmission_leak = 0);

        void initPriors();

        // Moves a sampled leak strictly inside its prior's range, as a new State starts it; a restored leak may come
        // from a Run with another prior.
        void keepLeakInsidePrior();

        // The Infection with `id`, and its position in `infections`. Throw std::out_of_range for an unknown id.
        [[nodiscard]] std::shared_ptr<InfectionEvent> infection(const std::string& id) const;
        [[nodiscard]] size_t indexOf(const std::string& id) const;

        // Visits every Quantity, in a fixed order: the Run scalars, then each Infection's values, each Locus's Allele
        // frequencies, each Infection's Latent genotype at each of its Loci (Loci in label order), then the Ordering.
        // The names are those of sampler targets, Trace columns and Checkpoint keys. A visitor handles every kind:
        //   scalar(name, parameter) · perInfection(name, index, parameter) · alleleFrequencies(locus, parameter)
        //   latentGenotype(index, locus, parameter) · ordering(ordering)
        template<typename Visitor>
        void quantities(Visitor&& v) const {
            v.scalar("mean_coi", meanCOI);
            v.scalar("mean_strains_tx", meanStrainsTransmitted);
            v.scalar("parent_set_size_prob", parentSetSizeProb);
            if (leakPrior) {
                v.scalar("transmission_leak", leakSize);
            }
            if (errorModel == ErrorModel::Hierarchical) {
                v.scalar("efp_mean", efpMean);
                v.scalar("efp_shape", efpShape);
                v.scalar("efn_mean", efnMean);
                v.scalar("efn_shape", efnShape);
            }
            for (size_t i = 0; i < infections.size(); ++i) {
                v.perInfection("expected_false_positives", i, expectedFalsePositives[i]);
                v.perInfection("expected_false_negatives", i, expectedFalseNegatives[i]);
                v.perInfection("infection_duration", i, infections[i]->infectionDuration());
            }
            for (const auto& [label, locus] : loci) {
                v.alleleFrequencies(locus, alleleFrequencies->alleleFrequencies(locus));
            }
            for (size_t i = 0; i < infections.size(); ++i) {
                for (const auto& [label, locus] : loci) {
                    if (infections[i]->latentGenotype().contains(locus)) {
                        v.latentGenotype(i, locus, infections[i]->latentGenotype(locus));
                    }
                }
            }
            v.ordering(infectionEventOrdering);
        }

        // The scalar and per-Infection Quantity names, in enumeration order.
        [[nodiscard]] std::vector<std::string> scalarNames() const;
        [[nodiscard]] std::vector<std::string> perInfectionNames() const;
        bool nullModel{};
        // The transmission leak ε as the Run recipe gave it: fixed, or with a leak prior the starting value.
        double transmissionLeak{};
        // The leak ε the Transmission memo reads: a Quantity (Run scalar transmission_leak) only with a leak prior.
        p_Parameterdouble leakSize;
        std::optional<LeakPrior> leakPrior{};
        ErrorModel errorModel = ErrorModel::Infection;
        std::map<std::string, std::shared_ptr<LocusImpl>> loci{};
        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        // Each infection's latent parent, summed out of the model. These name it in parent-set output (id `<id>_copy`)
        // and hold a fixed proposal anchor for the genotype samplers (see State.cpp); they are never sampled, read
        // by the likelihood, or checkpointed.
        std::vector<std::shared_ptr<InfectionEvent>> latentParents{};
        std::shared_ptr<core::containers::AllowedRelationships<InfectionEvent>> allowedRelationships;
        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        // Symptomatic vs Asymptomatic infection duration -- time between infection and detection

        p_DiscreteDist symptomaticInfectionDurationDist;
        p_DiscreteDist asymptomaticInfectionDurationDist;

        // Observation Process
        std::vector<p_Parameterdouble> expectedFalsePositives{};
        p_Parameterdouble obsFPRPriorShape;
        p_Parameterdouble obsFPRPriorScale;

        std::vector<p_Parameterdouble> expectedFalseNegatives{};
        p_Parameterdouble obsFNRPriorShape;
        p_Parameterdouble obsFNRPriorScale;

        // ErrorModel::Hierarchical: each rate's population mean and shape, sampled.
        p_Parameterdouble efpMean;
        p_Parameterdouble efpShape;
        p_Parameterdouble efnMean;
        p_Parameterdouble efnShape;

        // Parent Set Size Likelihood
        p_Parameterdouble parentSetSizeProb;
        p_Parameterdouble parentSetSizePriorAlpha;
        p_Parameterdouble parentSetSizePriorBeta;

        // Transmission term parameters
        p_Parameterdouble meanStrainsTransmitted;
        p_Parameterdouble meanStrainsTransmittedPriorShape;
        p_Parameterdouble meanStrainsTransmittedPriorScale;


        // Source Transmission Process
        p_Parameterdouble meanCOI;
        p_Parameterdouble meanCOIPriorShape;
        p_Parameterdouble meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::Model


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
