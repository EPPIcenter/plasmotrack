//
// Created by Maxwell Murphy on 10/29/21.
//

#include "State.h"
#include "core/io/serialize.h"

namespace transmission_nets::impl::ModelSix {
    State::State(const nlohmann::json& input) {
        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci);
        allowedParents = core::io::parseAllowedParentsFromJSON(input, infections);

        initPriors();

        alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        for (const auto& [locus_label, locus] : this->loci) {
            alleleFrequencies->addLocus(locus);
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

        for (size_t ii = 0; ii < infections.size(); ++ii) {
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(.1));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(.1));
        }

        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(.9);
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(.1);
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(5);
    }

    State::State(const nlohmann::json& input, const fs::path& outputDir) {
        // hotstart constructor
        auto paramOutputDir = outputDir / "parameters";
        auto epsPosFolder   = paramOutputDir / "eps_pos";
        auto epsNegFolder   = paramOutputDir / "eps_neg";
        auto infDurFolder   = paramOutputDir / "infection_duration";
        auto freqDir        = paramOutputDir / "allele_frequencies";
        auto genotypeDir    = paramOutputDir / "genotypes";

        loci           = core::io::parseLociFromJSON<LocusImpl>(input);
        infections     = core::io::parseInfectionsFromJSON<InfectionEvent, LocusImpl>(input, MAX_COI, loci);
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
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloaddouble(epsPosFolder / inf_file_name)));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(core::io::hotloaddouble(epsNegFolder / inf_file_name)));
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

        geometricGenerationProb = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "geo_gen_prob.csv"));
        lossProb                = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "loss_prob.csv"));
        meanCOI                 = std::make_shared<core::parameters::Parameter<double>>(core::io::hotloaddouble(paramOutputDir / "mean_coi.csv"));
    }

    void State::initPriors() {
        obsFPRPriorShape = std::make_shared<core::parameters::Parameter<double>>(1);
        obsFPRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.01);

        obsFNRPriorShape = std::make_shared<core::parameters::Parameter<double>>(1);
        obsFNRPriorScale = std::make_shared<core::parameters::Parameter<double>>(.01);

        geometricGenerationProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(1);
        geometricGenerationProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(1);

        lossProbPriorAlpha = std::make_shared<core::parameters::Parameter<double>>(10);
        lossProbPriorBeta  = std::make_shared<core::parameters::Parameter<double>>(90);

        meanCOIPriorShape = std::make_shared<core::parameters::Parameter<double>>(1);
        meanCOIPriorScale = std::make_shared<core::parameters::Parameter<double>>(5);

        infectionDurationShape = std::make_shared<core::parameters::Parameter<double>>(100);
        infectionDurationScale = std::make_shared<core::parameters::Parameter<double>>(1);
    }
}// namespace transmission_nets::impl::ModelSix
