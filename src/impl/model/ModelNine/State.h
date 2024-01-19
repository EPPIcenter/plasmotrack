//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_STATE_H
#define TRANSMISSION_NETWORKS_APP_STATE_H

#include "config.h"
#include "core/io/utils.h"
#include "core/io/parse_json.h"
#include "core/containers/AllowedRelationships.h"

#include <nlohmann/json.hpp>


namespace transmission_nets::impl::ModelNine {
    struct State {
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        using p_DiscreteDist = std::shared_ptr<core::distributions::DiscreteDistribution>;

        explicit State(const nlohmann::json& input,
                       const std::vector<core::computation::Probability>& symptomaticIDPDist,
                       const std::vector<core::computation::Probability>& asymptomaticIDPDist,
                       std::shared_ptr<boost::random::mt19937> rng,
                       bool null_model = false);
        State(const nlohmann::json& input,
              const std::vector<core::computation::Probability>& symptomaticIDPDist,
              const std::vector<core::computation::Probability>& asymptomaticIDPDist,
              std::shared_ptr<boost::random::mt19937> rng,
              const fs::path& outputDir, bool null_model = false);

        void initPriors();

        std::map<std::string, std::shared_ptr<LocusImpl>> loci{};
        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        std::vector<std::shared_ptr<InfectionEvent>> latentParents{};
        std::shared_ptr<core::containers::AllowedRelationships<InfectionEvent>> allowedRelationships;
        std::map<std::string, std::shared_ptr<ParentSetImpl>> parentSetList{};

        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        // Symptomatic vs Asymptomatic infection duration -- time between infection and detection
//        p_Parameterdouble symptomaticInfectionDurationShape;
//        p_Parameterdouble symptomaticInfectionDurationScale;
//        p_Parameterdouble asymptomaticInfectionDurationShape;
//        p_Parameterdouble asymptomaticInfectionDurationScale;

        p_DiscreteDist symptomaticInfectionDurationDist;
        p_DiscreteDist asymptomaticInfectionDurationDist;

        // Observation Process
        std::vector<p_Parameterdouble> expectedFalsePositives{};
        p_Parameterdouble obsFPRPriorShape;
        p_Parameterdouble obsFPRPriorScale;

        std::vector<p_Parameterdouble> expectedFalseNegatives{};
        p_Parameterdouble obsFNRPriorShape;
        p_Parameterdouble obsFNRPriorScale;

        // Node Transmission Process
        p_Parameterdouble meanStrainsTransmitted;
        p_Parameterdouble meanStrainsTransmittedPriorShape;
        p_Parameterdouble meanStrainsTransmittedPriorScale;
//        p_Parameterdouble geometricGenerationProb;
//        p_Parameterdouble geometricGenerationProbPriorAlpha;
//        p_Parameterdouble geometricGenerationProbPriorBeta;
//
//        p_Parameterdouble lossProb;
//        p_Parameterdouble lossProbPriorAlpha;
//        p_Parameterdouble lossProbPriorBeta;

//        p_Parameterdouble mutationProb;
//        p_Parameterdouble mutationProbPriorAlpha;
//        p_Parameterdouble mutationProbPriorBeta;

        // Source Transmission Process
        p_Parameterdouble meanCOI;
        p_Parameterdouble meanCOIPriorShape;
        p_Parameterdouble meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::ModelNine


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
