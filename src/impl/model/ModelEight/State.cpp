//
// Created by Maxwell Murphy on 10/29/21.
//

#include "State.h"
#include "core/io/serialize.h"

namespace transmission_nets::impl::ModelEight {
    State::State(const nlohmann::json& input, std::shared_ptr<boost::random::mt19937> rng) {
        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, rng);
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
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(.1);
        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(.05);
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(1.01);
    }

    State::State(const nlohmann::json& input, std::shared_ptr<boost::random::mt19937> rng, const fs::path& outputDir) {
        // hotstart constructor
        auto paramOutputDir = outputDir / "parameters";
        auto epsPosFolder   = paramOutputDir / "eps_pos";
        auto epsNegFolder   = paramOutputDir / "eps_neg";
        auto infDurFolder   = paramOutputDir / "infection_duration";
        auto freqDir        = paramOutputDir / "allele_frequencies";
        auto genotypeDir    = paramOutputDir / "genotypes";
        auto latentParentsDir = paramOutputDir / "latent_parents";

        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci, rng);
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
            infection->infectionDuration()->initializeValue(core::io::hotloadDouble(infDurFolder / inf_file_name));
            for (const auto& [label, locus] : loci) {
                infection->latentGenotype(locus)->initializeValue(GeneticsImpl(core::io::hotloadString(infDir / (label + ".csv"))));
            }

            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloadDouble(epsPosFolder / inf_file_name)));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloadDouble(epsNegFolder / inf_file_name)));
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

        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "geo_gen_prob.csv"));
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "loss_prob.csv"));
        mutationProb            = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "mutation_prob.csv"));
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloadDouble(paramOutputDir / "mean_coi.csv"));
    }

    void State::initPriors() {
        obsFPRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFPRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        obsFNRPriorShape = std::make_shared<core::parameters::Parameter<double>>(10);
        obsFNRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.001);

        geometricGenerationProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(10);
        geometricGenerationProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(1);

        lossProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
        lossProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(1);

//        mutationProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
//        mutationProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(99999);

        meanCOIPriorShape = std::make_shared<core::parameters::Parameter<double>>(20);
        meanCOIPriorScale = std::make_shared<core::parameters::Parameter<double>>(.1);

        infectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(10);
        infectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(10);

    }
}// namespace transmission_nets::impl::ModelEight
