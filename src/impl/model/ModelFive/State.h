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


namespace transmission_nets::impl::ModelFive {
    struct State {
        //        State(std::map<std::string, LocusImpl *> loci, const std::vector<InfectionEvent *> &infections, std::map<InfectionEvent *, std::vector<InfectionEvent *>> allowedParents);

        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        explicit State(const nlohmann::json& input);
        State(const nlohmann::json& input, const fs::path& outputDir);

        std::map<std::string, std::shared_ptr<LocusImpl>> loci{};
        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents{};

        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        p_Parameterdouble infectionDurationShape;
        p_Parameterdouble infectionDurationScale;


        // Observation Process
        std::vector<p_Parameterdouble> observationFalsePositiveRates{};
        p_Parameterdouble obsFPRPriorAlpha;
        p_Parameterdouble obsFPRPriorBeta;

        std::vector<p_Parameterdouble> observationFalseNegativeRates{};
        p_Parameterdouble obsFNRPriorAlpha;
        p_Parameterdouble obsFNRPriorBeta;

        // Node Transmission Process
        p_Parameterdouble geometricGenerationProb;
        p_Parameterdouble geometricGenerationProbPriorAlpha;
        p_Parameterdouble geometricGenerationProbPriorBeta;

        p_Parameterdouble lossProb;
        p_Parameterdouble lossProbPriorAlpha;
        p_Parameterdouble lossProbPriorBeta;


        p_Parameterdouble mutationProb;
        p_Parameterdouble mutationProbPriorAlpha;
        p_Parameterdouble mutationProbPriorBeta;

        // Source Transmission Process
        p_Parameterdouble meanCOI;
        p_Parameterdouble meanCOIPriorShape;
        p_Parameterdouble meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::ModelFive


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
