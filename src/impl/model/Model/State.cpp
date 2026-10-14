
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

    }
    namespace {
        // Records the names the enumeration gives, once each.
        struct Names {
            std::vector<std::string> scalars;
            std::vector<std::string> perInfectionQuantities;

            template<typename P>
            void scalar(const char* name, const P&) { scalars.emplace_back(name); }
            template<typename P>
            void perInfection(const char* name, const size_t index, const P&) {
                if (index == 0) {
                    perInfectionQuantities.emplace_back(name);
                }
            }
            template<typename L, typename P>
            void alleleFrequencies(const L&, const P&) {}
            template<typename L, typename P>
            void latentGenotype(size_t, const L&, const P&) {}
            template<typename O>
            void ordering(const O&) {}
        };
    }// namespace

    std::vector<std::string> State::scalarNames() const {
        Names names;
        quantities(names);
        return names.scalars;
    }

    std::vector<std::string> State::perInfectionNames() const {
        Names names;
        quantities(names);
        return names.perInfectionQuantities;
    }
}// namespace transmission_nets::impl::Model
