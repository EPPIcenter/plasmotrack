//
// Created by Maxwell Murphy on 10/29/21.
//

#include "State.h"
#include "core/io/serialize.h"

namespace transmission_nets::impl::Model {
    State::State(
            const nlohmann::json& input,
            const std::vector<core::computation::Probability>& symptomaticIDPDist,
            const std::vector<core::computation::Probability>& asymptomaticIDPDist,
            const std::shared_ptr<boost::random::mt19937>& rng,
            const bool null_model) {
        null_model_ = null_model;
        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, rng, null_model);
        allowedRelationships = core::io::parseAllowedParentsFromJSON(input, infections);

        initPriors();

        // alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        // for (const auto& [locus_label, locus] : this->loci) {
        //     alleleFrequencies->addLocus(locus);
        // }

        alleleFrequencies = core::io::parseAlleleFrequenciesFromJSON<AlleleFrequencyContainerImpl>(input, loci, .05);

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

        for (const auto& infection : infections) {
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(.01));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(.01));
            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
        }

//        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(.95);
//        lossProb                = std::make_shared<core::parameters::Parameter<double>>(.5);
//        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(.05);
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(1.01);
        meanStrainsTransmitted  = std::make_shared<core::parameters::Parameter<double>>(2.00);
        parentSetSizeProb = std::make_shared<core::parameters::Parameter<double>>(.9);
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
        null_model_ = null_model;
        auto paramOutputDir = outputDir / "parameters";
        auto epsPosFolder   = paramOutputDir / "eps_pos";
        auto epsNegFolder   = paramOutputDir / "eps_neg";
        auto infDurFolder   = paramOutputDir / "infection_duration";
        auto freqDir        = paramOutputDir / "allele_frequencies";
        auto genotypeDir    = paramOutputDir / "genotypes";
        auto latentParentsDir = paramOutputDir / "latent_parents";

        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, std::move(rng), null_model);
        allowedRelationships = core::io::parseAllowedParentsFromJSON(input, infections);

        initPriors();

        // alleleFrequencies = core::io::parseAlleleFrequenciesFromJSON<AlleleFrequencyContainerImpl>(input, loci);
        alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        for (const auto& locus : this->loci | std::views::values) {
            alleleFrequencies->addLocus(locus);
            auto hotloadFreq = core::datatypes::Simplex(core::io::hotloadVector(freqDir / (locus->label + ".csv.gz")));
            alleleFrequencies->alleleFrequencies(locus)->initializeValue(std::move(hotloadFreq));
        }

        for (auto& infection : infections) {
            auto infDir        = genotypeDir / core::io::makePathValid(infection->id());
            auto inf_file_name = core::io::makePathValid(infection->id()) + ".csv.gz";
            infection->infectionDuration()->initializeValue(core::io::hotloadDouble(infDurFolder / inf_file_name));
            for (const auto& [label, locus] : loci) {
                infection->latentGenotype(locus)->initializeValue(GeneticsImpl(core::io::hotloadString(infDir / (label + ".csv.gz"))));
            }

            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloadDouble(epsPosFolder / inf_file_name)));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloadDouble(epsNegFolder / inf_file_name)));
        }

        for (auto& infection : latentParents) {
            auto infDir        = latentParentsDir / core::io::makePathValid(infection->id());
            for (const auto& [label, locus] : loci) {
                infection->latentGenotype(locus)->initializeValue(GeneticsImpl(core::io::hotloadString(infDir / (label + ".csv.gz"))));
            }
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

//        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "geo_gen_prob.csv"));
//        lossProb                = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "loss_prob.csv"));
//        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "mutation_prob.csv"));

        meanCOI = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "mean_coi.csv.gz"));
        meanStrainsTransmitted = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "mean_strains_tx.csv.gz"));
        parentSetSizeProb = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "parent_set_size_prob.csv.gz"));
        symptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(symptomaticIDPDist);
        asymptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(asymptomaticIDPDist);
    }

    void State::initPriors() {
        obsFPRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFPRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        obsFNRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFNRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        meanStrainsTransmittedPriorShape = std::make_shared<core::parameters::Parameter<double>>(20);
        meanStrainsTransmittedPriorScale = std::make_shared<core::parameters::Parameter<double>>(.1);

        meanCOIPriorShape = std::make_shared<core::parameters::Parameter<double>>(20);
        meanCOIPriorScale = std::make_shared<core::parameters::Parameter<double>>(.1);

        parentSetSizePriorAlpha = std::make_shared<core::parameters::Parameter<double>>(10);
        parentSetSizePriorBeta = std::make_shared<core::parameters::Parameter<double>>(1);

        /*
         * Symptomatic vs Asymptomatic Infection Duration
         * Parameterizes the time between infection and clinical detection
         */
//        symptomaticInfectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(10);
//        symptomaticInfectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(10);
//        asymptomaticInfectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(1);
//        asymptomaticInfectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(100);
    }
}// namespace transmission_nets::impl::ModelNine
