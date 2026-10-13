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

        alleleFrequencies = core::io::parseAlleleFrequenciesFromJSON<AlleleFrequencyContainerImpl>(input, loci, MIN_ALLELE_FREQUENCY);

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(this->infections);

        for (const auto& infection : infections) {
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(.01));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(.01));
            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
            // The anchor: at each locus the observed genotype, or every allele where nothing was observed. A function
            // of the data alone, so a resumed Run proposes as the original would.
            for (const auto& locus : latentParents.back()->loci()) {
                GeneticsImpl anchor(static_cast<int>(locus->totalAlleles()));
                anchor.reset();
                if (infection->observedGenotype().contains(locus) and infection->observedGenotype(locus)->value().totalPositiveCount() > 0) {
                    anchor = infection->observedGenotype(locus)->value();
                } else {
                    for (unsigned int a = 0; a < locus->totalAlleles(); ++a) {
                        anchor.set(a);
                    }
                }
                latentParents.back()->latentGenotype(locus)->initializeValue(anchor);
            }
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
