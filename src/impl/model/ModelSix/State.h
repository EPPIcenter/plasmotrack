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


namespace transmission_nets::impl::ModelSix {
    struct State {
        //        State(std::map<std::string, LocusImpl *> loci, const std::vector<InfectionEvent *> &infections, std::map<InfectionEvent *, std::vector<InfectionEvent *>> allowedParents);

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
        explicit State(const nlohmann::json &input);
        State(const nlohmann::json &input, const fs::path &outputDir);

        std::map<std::string, std::shared_ptr<LocusImpl>> loci{};
        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents{};

        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        p_ParameterDouble infectionDurationShape;
        p_ParameterDouble infectionDurationScale;

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
}// namespace transmission_nets::impl::ModelSix


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
