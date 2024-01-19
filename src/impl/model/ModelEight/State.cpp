//
// Created by Maxwell Murphy on 10/29/21.
//

#include "State.h"
#include "core/io/serialize.h"

namespace transmission_nets::impl::ModelEight {
    State::State(
            const nlohmann::json& input,
            const std::vector<core::computation::Probability>& symptomaticIDPDist,
            const std::vector<core::computation::Probability>& asymptomaticIDPDist,
            std::shared_ptr<boost::random::mt19937> rng,
            const bool null_model) {
        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, rng, null_model);
        allowedParents = core::io::parseAllowedParentsFromJSON(input, infections);

        initPriors();

        alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        for (const auto& [locus_label, locus] : this->loci) {
            alleleFrequencies->addLocus(locus);
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);


        for (const auto& infection : infections) {
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(.01));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(.01));
//            latentParents.push_back(std::make_shared<InfectionEvent>(*infection, "", false));
            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
        }

        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(.95);
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(.5);
        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(.05);
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(1.01);
        symptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(symptomaticIDPDist);
        asymptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(asymptomaticIDPDist);
    }

    State::State(
            const nlohmann::json& input,
            const std::vector<core::computation::Probability>& symptomaticIDPDist,
            const std::vector<core::computation::Probability>& asymptomaticIDPDist,
            std::shared_ptr<boost::random::mt19937> rng,
            const fs::path& outputDir,
            const bool null_model) {
        // hotstart constructor
        auto paramOutputDir = outputDir / "parameters";
        auto epsPosFolder   = paramOutputDir / "eps_pos";
        auto epsNegFolder   = paramOutputDir / "eps_neg";
        auto infDurFolder   = paramOutputDir / "infection_duration";
        auto freqDir        = paramOutputDir / "allele_frequencies";
        auto genotypeDir    = paramOutputDir / "genotypes";
        auto latentParentsDir = paramOutputDir / "latent_parents";

        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, rng, null_model);
        allowedParents = core::io::parseAllowedParentsFromJSON(input, infections);

        initPriors();

        alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        for (const auto& [locus_label, locus] : this->loci) {
            alleleFrequencies->addLocus(locus);
            auto hotloadFreq = core::datatypes::Simplex(core::io::hotloadVector(freqDir / (locus->label + ".csv")));
            alleleFrequencies->alleleFrequencies(locus)->initializeValue(hotloadFreq);
        }

        for (auto& infection : infections) {
            auto infDir        = genotypeDir / core::io::makePathValid(infection->id());
            auto inf_file_name = core::io::makePathValid(infection->id()) + ".csv";
            infection->infectionDuration()->initializeValue(core::io::hotloaddouble(infDurFolder / inf_file_name));
            for (const auto& [label, locus] : loci) {
                infection->latentGenotype(locus)->initializeValue(GeneticsImpl(core::io::hotloadString(infDir / (label + ".csv"))));
            }

            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloaddouble(epsPosFolder / inf_file_name)));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloaddouble(epsNegFolder / inf_file_name)));
        }

        for (auto& infection : latentParents) {
            auto infDir        = latentParentsDir / core::io::makePathValid(infection->id());
            auto inf_file_name = core::io::makePathValid(infection->id()) + ".csv";
            for (const auto& [label, locus] : loci) {
                infection->latentGenotype(locus)->initializeValue(GeneticsImpl(core::io::hotloadString(infDir / (label + ".csv"))));
            }
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "geo_gen_prob.csv"));
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "loss_prob.csv"));
        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "mutation_prob.csv"));
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "mean_coi.csv"));
        symptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(symptomaticIDPDist);
        asymptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(asymptomaticIDPDist);
    }

    void State::initPriors() {
        obsFPRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFPRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        obsFNRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFNRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        geometricGenerationProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
        geometricGenerationProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(1);

        lossProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
        lossProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(1);

//        mutationProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
//        mutationProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(99999);

        meanCOIPriorShape = std::make_shared<core::parameters::Parameter<double>>(20);
        meanCOIPriorScale = std::make_shared<core::parameters::Parameter<double>>(.1);

        /*
         * Symptomatic vs Asymptomatic Infection Duration
         * Parameterizes the time between infection and clinical detection
         */
//        symptomaticInfectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(10);
//        symptomaticInfectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(10);
//        asymptomaticInfectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(1);
//        asymptomaticInfectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(100);

    }
}// namespace transmission_nets::impl::ModelEight
