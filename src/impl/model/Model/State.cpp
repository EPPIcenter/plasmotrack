//
// The State's construction from a Study: the data rules and the refusals (see State.h).
//

#include "State.h"

#include "core/io/serialize.h"

#include <fmt/core.h>

#include <algorithm>
#include <stdexcept>

namespace transmission_nets::impl::Model {
    namespace {
        [[noreturn]] void refuseInput(std::string message) {
            throw Refused(refuse(RefusalKind::InvalidInput, std::move(message)));
        }

        // An all-zero (or empty) genotype is unobserved.
        bool missingGenotype(const std::string& alleles) {
            return std::ranges::all_of(alleles, [](const char s) { return s == '0'; });
        }

        // The genotype with fewer than `max_coi` Alleles present: the first `max_coi - 1`, in Allele order. The cap
        // is one below the limit, as it has always been; the Golden runs depend on it.
        std::string constrainCOI(const std::string& genotype, const int max_coi) {
            std::string out;
            int total_coi = 0;
            for (const auto& allele : genotype) {
                if (allele == '1') {
                    total_coi += 1;
                }
                if (total_coi < max_coi) {
                    out += allele;
                } else {
                    out += '0';
                }
            }
            return out;
        }
    }// namespace

    State::State(const Study& study, const std::shared_ptr<boost::random::mt19937>& rng, const bool null_model) {
        nullModel = null_model;

        // Loci, in input order, so that their uids follow the input.
        std::vector<std::shared_ptr<LocusImpl>> lociInOrder;
        for (const auto& locus : study.loci) {
            if (locus.alleleFrequencies.size() != locus.alleles) {
                refuseInput(fmt::format("Locus {} has {} Allele frequencies for {} Alleles.", locus.label, locus.alleleFrequencies.size(), locus.alleles));
            }
            lociInOrder.push_back(std::make_shared<LocusImpl>(locus.label, static_cast<int>(locus.alleles)));
            if (!loci.emplace(locus.label, lociInOrder.back()).second) {
                refuseInput(fmt::format("Locus {} is listed twice.", locus.label));
            }
        }

        // Infections, in input order. An unobserved Locus has no Observed genotype, and its Latent genotype starts as
        // one random Allele; an observed one starts as the observation, capped at MAX_COI Alleles. Under the null
        // model every Locus is unobserved and every observation time the same.
        for (const auto& described : study.infections) {
            if (described.observationTime < 0) {
                refuseInput(fmt::format("Observation time for Infection {} is negative.", described.id));
            }
            const double observationTime = null_model ? 1000 : described.observationTime;
            auto infection               = std::make_shared<InfectionEvent>(described.id, observationTime, described.symptomatic);
            auto addUnobserved           = [&](const std::shared_ptr<LocusImpl>& locus) {
                std::string latent_genetics;
                latent_genetics.resize(locus->totalAlleles(), '0');
                boost::random::uniform_int_distribution<> dist(0, static_cast<int>(locus->totalAlleles()) - 1);
                latent_genetics[dist(*rng)] = '1';
                infection->addLatentGenetics(locus, latent_genetics);
            };
            for (const auto& genotype : described.observedGenotypes) {
                const auto locus = loci.find(genotype.locus);
                if (locus == loci.end()) {
                    refuseInput(fmt::format("Locus {} for Infection {} does not exist.", genotype.locus, described.id));
                }
                if (genotype.alleles.size() != locus->second->totalAlleles()) {
                    refuseInput(fmt::format("Infection {} has a genotype of {} Alleles at Locus {}, which has {}.", described.id, genotype.alleles.size(), genotype.locus, locus->second->totalAlleles()));
                }
                if (!null_model and !missingGenotype(genotype.alleles)) {
                    infection->addGenetics(locus->second, genotype.alleles, constrainCOI(genotype.alleles, MAX_COI));
                } else {
                    addUnobserved(locus->second);
                }
            }
            // A Locus missing from the list is unobserved too. The model assumes every Infection has every Locus:
            // parents' and children's genotypes are read at each other's Loci.
            for (const auto& [label, locus] : loci) {
                if (not infection->latentGenotype().contains(locus)) {
                    addUnobserved(locus);
                }
            }
            infections.push_back(std::move(infection));
        }

        // The parents the study permits each Infection, which with the Ordering give its Candidate parents.
        allowedRelationships = std::make_shared<core::containers::AllowedRelationships<InfectionEvent>>();
        for (size_t i = 0; i < infections.size(); ++i) {
            allowedRelationships->addInfectionEvent(infections[i]);
            for (const auto& parentId : study.infections[i].permittedParents) {
                const auto parent = std::ranges::find_if(infections, [&](const auto& candidate) { return candidate->id() == parentId; });
                if (parent == infections.end()) {
                    refuseInput(fmt::format("Infection {} has a permitted parent {} that does not exist.", infections[i]->id(), parentId));
                }
                allowedRelationships->addParent(infections[i], *parent);
                allowedRelationships->addChild(*parent, infections[i]);
            }
        }

        initPriors();

        // Allele frequencies, floored at MIN_ALLELE_FREQUENCY.
        alleleFrequencies = std::make_shared<AlleleFrequencyContainerImpl>();
        for (size_t l = 0; l < study.loci.size(); ++l) {
            alleleFrequencies->addLocus(lociInOrder[l]);
            auto frequencies = study.loci[l].alleleFrequencies;
            for (auto& frequency : frequencies) {
                if (frequency < MIN_ALLELE_FREQUENCY) {
                    frequency = MIN_ALLELE_FREQUENCY;
                }
            }
            alleleFrequencies->alleleFrequencies(lociInOrder[l])->initializeValue(AlleleFrequencyImpl(frequencies));
        }

        infectionEventOrdering = std::make_shared<OrderingImpl>();
        infectionEventOrdering->addElements(infections, [&](const auto& infection) { return allowedRelationships->permittedParents(infection); });

        for (const auto& infection : infections) {
            expectedFalsePositives.emplace_back(new core::parameters::Parameter<double>(.01));
            expectedFalseNegatives.emplace_back(new core::parameters::Parameter<double>(.01));
            latentParents.push_back(std::make_shared<InfectionEvent>(*infection));
            // The anchor: at each locus the observed genotype, or every allele where nothing was observed. A function
            // of the data alone, so a resumed Run proposes as the original would.
            for (const auto& locus : latentParents.back()->loci()) {
                Genotype anchor(static_cast<int>(locus->totalAlleles()));
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

        meanCOI                           = std::make_shared<core::parameters::Parameter<double>>(1.01);
        meanStrainsTransmitted            = std::make_shared<core::parameters::Parameter<double>>(2.00);
        parentSetSizeProb                 = std::make_shared<core::parameters::Parameter<double>>(.9);
        symptomaticInfectionDurationDist  = std::make_shared<core::distributions::DiscreteDistribution>(study.symptomaticInfectionDurations);
        asymptomaticInfectionDurationDist = std::make_shared<core::distributions::DiscreteDistribution>(study.asymptomaticInfectionDurations);
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

    size_t State::indexOf(const std::string& id) const {
        for (size_t i = 0; i < infections.size(); ++i) {
            if (infections[i]->id() == id) {
                return i;
            }
        }
        throw std::out_of_range("no Infection " + id);
    }

    std::shared_ptr<InfectionEvent> State::infection(const std::string& id) const {
        return infections[indexOf(id)];
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
