
#ifndef TRANSMISSION_NETWORKS_APP_STATE_H
#define TRANSMISSION_NETWORKS_APP_STATE_H

#include "config.h"
#include "core/io/utils.h"
#include "core/io/parse_json.h"
#include "core/containers/AllowedRelationships.h"

#include <nlohmann/json.hpp>


namespace transmission_nets::impl::Model {
    struct State {
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        using p_DiscreteDist = std::shared_ptr<core::distributions::DiscreteDistribution>;

        explicit State(const nlohmann::json& input,
                       const std::vector<core::computation::Probability>& symptomaticIDPDist,
                       const std::vector<core::computation::Probability>& asymptomaticIDPDist,
                       const std::shared_ptr<boost::random::mt19937>& rng,
                       bool null_model = false);

        void initPriors();

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
        bool null_model_{};
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

        // Parent Set Size Likelihood
        p_Parameterdouble parentSetSizeProb;
        p_Parameterdouble parentSetSizePriorAlpha;
        p_Parameterdouble parentSetSizePriorBeta;

        // Node Transmission Process
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
