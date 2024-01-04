//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_STATE_H
#define TRANSMISSION_NETWORKS_APP_STATE_H

#include "config.h"
#include "core/datatypes/Simplex.h"
#include "core/io/parse_json.h"
#include "core/io/utils.h"

#include <nlohmann/json.hpp>
#include <random>


namespace transmission_nets::impl::ModelEight {
    struct State {
        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
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
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents{};
        std::map<std::string, std::shared_ptr<ParentSetImpl>> parentSetList{};

        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        // Symptomatic vs Asymptomatic infection duration -- time between infection and detection
//        p_ParameterDouble symptomaticInfectionDurationShape;
//        p_ParameterDouble symptomaticInfectionDurationScale;
//        p_ParameterDouble asymptomaticInfectionDurationShape;
//        p_ParameterDouble asymptomaticInfectionDurationScale;

        p_DiscreteDist symptomaticInfectionDurationDist;
        p_DiscreteDist asymptomaticInfectionDurationDist;

        // Observation Process
        std::vector<p_ParameterDouble> expectedFalsePositives{};
        p_ParameterDouble obsFPRPriorShape;
        p_ParameterDouble obsFPRPriorScale;

        std::vector<p_ParameterDouble> expectedFalseNegatives{};
        p_ParameterDouble obsFNRPriorShape;
        p_ParameterDouble obsFNRPriorScale;

        // Node Transmission Process
        p_ParameterDouble geometricGenerationProb;
        p_ParameterDouble geometricGenerationProbPriorAlpha;
        p_ParameterDouble geometricGenerationProbPriorBeta;

        p_ParameterDouble lossProb;
        p_ParameterDouble lossProbPriorAlpha;
        p_ParameterDouble lossProbPriorBeta;

        p_ParameterDouble mutationProb;
        p_ParameterDouble mutationProbPriorAlpha;
        p_ParameterDouble mutationProbPriorBeta;

        // Source Transmission Process
        p_ParameterDouble meanCOI;
        p_ParameterDouble meanCOIPriorShape;
        p_ParameterDouble meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::ModelEight


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
