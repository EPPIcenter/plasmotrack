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

        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;
        explicit State(const nlohmann::json& input);
        State(const nlohmann::json& input, const fs::path& outputDir);

        void initPriors();

        std::map<std::string, std::shared_ptr<LocusImpl>> loci{};
        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents{};

        std::shared_ptr<AlleleFrequencyContainerImpl> alleleFrequencies;

        // Network Structure
        std::shared_ptr<OrderingImpl> infectionEventOrdering;

        p_Parameterfloat infectionDurationShape;
        p_Parameterfloat infectionDurationScale;

        // Observation Process
        std::vector<p_Parameterfloat> expectedFalsePositives{};
        p_Parameterfloat obsFPRPriorShape;
        p_Parameterfloat obsFPRPriorScale;

        std::vector<p_Parameterfloat> expectedFalseNegatives{};
        p_Parameterfloat obsFNRPriorShape;
        p_Parameterfloat obsFNRPriorScale;

        // Node Transmission Process
        p_Parameterfloat geometricGenerationProb;
        p_Parameterfloat geometricGenerationProbPriorAlpha;
        p_Parameterfloat geometricGenerationProbPriorBeta;

        p_Parameterfloat lossProb;
        p_Parameterfloat lossProbPriorAlpha;
        p_Parameterfloat lossProbPriorBeta;

//        p_Parameterfloat mutationProb;
//        p_Parameterfloat mutationProbPriorAlpha;
//        p_Parameterfloat mutationProbPriorBeta;

        // Source Transmission Process
        p_Parameterfloat meanCOI;
        p_Parameterfloat meanCOIPriorShape;
        p_Parameterfloat meanCOIPriorScale;
    };
}// namespace transmission_nets::impl::ModelSix


#endif//TRANSMISSION_NETWORKS_APP_STATE_H
