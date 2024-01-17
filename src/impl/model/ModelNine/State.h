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
        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;
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
//        p_Parameterfloat symptomaticInfectionDurationShape;
//        p_Parameterfloat symptomaticInfectionDurationScale;
//        p_Parameterfloat asymptomaticInfectionDurationShape;
//        p_Parameterfloat asymptomaticInfectionDurationScale;

        p_DiscreteDist symptomaticInfectionDurationDist;
        p_DiscreteDist asymptomaticInfectionDurationDist;

        // Observation Process
        std::vector<p_Parameterfloat> expectedFalsePositives{};
        p_Parameterfloat obsFPRPriorShape;
        p_Parameterfloat obsFPRPriorScale;

        std::vector<p_Parameterfloat> expectedFalseNegatives{};
        p_Parameterfloat obsFNRPriorShape;
        p_Parameterfloat obsFNRPriorScale;

        // Node Transmission Process
        p_Parameterfloat meanStrainsTransmitted;
        p_Parameterfloat meanStrainsTransmittedPriorShape;
        p_Parameterfloat meanStrainsTransmittedPriorScale;
//        p_Parameterfloat geometricGenerationProb;
//        p_Parameterfloat geometricGenerationProbPriorAlpha;
//        p_Parameterfloat geometricGenerationProbPriorBeta;
//
//        p_Parameterfloat lossProb;
//        p_Parameterfloat lossProbPriorAlpha;
//        p_Parameterfloat lossProbPriorBeta;

//        p_Parameterfloat mutationProb;
//        p_Parameterfloat mutationProbPriorAlpha;
//        p_Parameterfloat mutationProbPriorBeta;

        // Source Transmission Process
        p_Parameterfloat meanCOI;
        p_Parameterfloat meanCOIPriorShape;
        p_Parameterfloat meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::ModelNine


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
