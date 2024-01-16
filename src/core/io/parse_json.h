//
// Created by Maxwell Murphy on 5/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARSE_JSON_H
#define TRANSMISSION_NETWORKS_APP_PARSE_JSON_H

#include "core/containers/AllowedRelationships.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/datatypes/Simplex.h"
#include "core/distributions/DiscreteDistribution.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <boost/type.hpp>
#include <iostream>
#include <memory>

namespace transmission_nets::core::io {
    using nlohmann::json;

    inline json loadJSON(std::istream& input) {
        fmt::print("Loading JSON");
        auto j = json::parse(input);
        fmt::print("...done\n");
        return j;
    }

    inline bool missingGenotype(const std::string& alleleStr) {
        bool isEmpty = alleleStr.empty();
        bool isMissing = std::all_of(alleleStr.begin(), alleleStr.end(), [](const char s) { return s == '0'; });
        return isEmpty or isMissing;
    }

    inline std::string constrainCOI(const std::string& genotype, int max_coi) {
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


    template<typename LocusImpl>
    std::map<std::string, std::shared_ptr<LocusImpl>> parseLociFromJSON(
            const json& input,
            const char lociKey[] = "loci",
            const char locusLabelKey[] = "locus",
            const char numAllelesKey[] = "num_alleles") {

        std::map<std::string, std::shared_ptr<LocusImpl>> locusMap{};

        for (const auto& loc : input.at(lociKey)) {
            const std::string locus_label = loc.at(locusLabelKey);
            int num_alleles = loc.at(numAllelesKey);
            locusMap.emplace(locus_label, std::make_shared<containers::Locus>(locus_label, num_alleles));
        }

        return locusMap;
    }

    template<typename AlleleFrequencyContainerImpl, typename LocusImpl>
    std::shared_ptr<AlleleFrequencyContainerImpl> parseAlleleFrequenciesFromJSON(
            const json& input,
            std::map<std::string, std::shared_ptr<LocusImpl>> loci,
            const double minAlleleFreq = .01,
            const char alleleFrequencyKey[] = "allele_freqs",
            const char lociKey[] = "loci",
            const char locusLabelKey[] = "locus") {
        auto afContainer = std::make_shared<AlleleFrequencyContainerImpl>();

        for (const auto& afEntry : input.at(lociKey)) {
            std::string locusLabel = afEntry.at(locusLabelKey);
            auto locus = loci.at(locusLabel);
            afContainer->addLocus(locus);
            std::vector<double> afToLoad = afEntry.at(alleleFrequencyKey);

            for (auto& af : afToLoad) {
                if (af < minAlleleFreq) {
                    af = minAlleleFreq;
                }
            }

            auto af = core::datatypes::Simplex(std::vector<double>(afToLoad));
            afContainer->alleleFrequencies(locus)->initializeValue(af);
        }

        return afContainer;
    }

    template<typename InfectionEvent, typename LocusImpl, typename Engine = boost::random::mt19937>
    std::vector<std::shared_ptr<InfectionEvent>> parseInfectionsFromJSON(
            const json& input,
            const int max_coi,
            std::map<std::string, std::shared_ptr<LocusImpl>> loci,
            std::shared_ptr<Engine> rng,
            const bool null_model = false,
            const char infectionsKey[] = "nodes",
            const char obsGenotypesKey[] = "observed_genotype",
            const char idKey[] = "id",
            const char observationDateKey[] = "observation_time",
            const char symptomaticKey[] = "symptomatic",
            const char genotypeKey[] = "genotype",
            const char locusKey[] = "locus") {

        if (null_model) {
            std::cout << "WARNING: Null model is enabled. All observed genotypes and infection times will be ignored.\n";
        }

        std::vector<std::shared_ptr<InfectionEvent>> infections{};
        for (const auto& inf : input.at(infectionsKey)) {
            bool symptomatic = true;
            if (inf.count(symptomaticKey) != 0) {
                symptomatic = inf.at(symptomaticKey);
            }

            double obs_time = inf.at(observationDateKey);
            if (obs_time < 0) {
                std::cerr << "Observation time for node " << inf.at(idKey) << " is negative.\n";
                exit(1);
            }

            if (null_model) {
                obs_time = 1000;
            }

            infections.push_back(std::make_shared<InfectionEvent>(inf.at(idKey), obs_time, symptomatic));
            for (const auto& genetics : inf.at(obsGenotypesKey)) {
                auto locusLabel = genetics.at(locusKey);
                auto locusItr = loci.find(locusLabel);
                if (locusItr == loci.end()) {
                    std::cerr << "Locus " << locusLabel << " for node " << infections.back()->id() << " does not exist.\n";
                    exit(1);
                }
                if (!null_model and !missingGenotype(genetics.at(genotypeKey))) {
                    std::string obs_genetics = genetics.at(genotypeKey);
                    infections.back()->addGenetics(locusItr->second, obs_genetics, constrainCOI(obs_genetics, max_coi));
                } else {
                    std::string latent_genetics;
                    latent_genetics.resize(locusItr->second->totalAlleles(), '0');
                    boost::random::uniform_int_distribution<> dist(0, locusItr->second->totalAlleles() - 1);
                    latent_genetics[dist(*rng)] = '1';
                    infections.back()->addLatentGenetics(locusItr->second, latent_genetics);
                }
            }
        }

        return infections;
    }


    template<typename InfectionEvent>
    containers::AllowedRelationships<InfectionEvent> parseAllowedParentsFromJSON(
            const json& input,
            std::vector<std::shared_ptr<InfectionEvent>> infections,
            const char infectionsKey[] = "nodes",
            const char idKey[] = "id",
            const char allowedParentsKey[] = "allowed_parents") {

        // std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents{};
        // std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedChildren{};

        containers::AllowedRelationships<InfectionEvent> allowedRelationships{};

        for (const auto& targetInfection : infections) {
            for (const auto& inf : input.at(infectionsKey)) {
                auto infectionId = inf.at(idKey);
                if (infectionId == targetInfection->id()) {
                    for (const auto& parentId : inf.at(allowedParentsKey)) {
                        auto parentInf = std::ranges::find_if(
                                infections,
                                [&parentId](const std::shared_ptr<InfectionEvent> candidateInf) {
                                    return candidateInf->id() == parentId;
                                });

                        allowedRelationships.allowedParents[targetInfection].push_back(*parentInf);
                        allowedRelationships.allowedChildren[*parentInf].push_back(targetInfection);
                    }
                    break;
                }
            }
        }

        return allowedRelationships;
    }
}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_PARSE_JSON_H
