//
// Created by Maxwell Murphy on 5/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARSE_JSON_H
#define TRANSMISSION_NETWORKS_APP_PARSE_JSON_H

#include <nlohmann/json.hpp>
#include <iostream>

#include "core/containers/Locus.h"
#include "core/containers/Infection.h"

namespace transmission_nets::core::io {
    using nlohmann::json;

    inline json loadJSON(std::istream &input) {
        auto j = json::parse(input);
        return j;
    }

    inline bool missingGenotype(const std::string& alleleStr) {
        bool isEmpty = alleleStr.empty();
        bool isMissing = std::all_of(alleleStr.begin(), alleleStr.end(), [](const char s) { return s == '0'; });
        return isEmpty or isMissing;
    }


    template <typename LocusImpl>
    std::map<std::string, LocusImpl*> parseLociFromJSON(
            const json &input,
            const char lociKey[] = "loci",
            const char locusLabelKey[] = "locus",
            const char numAllelesKey[] = "num_alleles") {

        std::map<std::string, LocusImpl*> locusMap{};

        for(const auto& loc : input.at(lociKey)) {
            const std::string locus_label = loc.at(locusLabelKey);
            int num_alleles = loc.at(numAllelesKey);
            locusMap.emplace(locus_label, new containers::Locus(locus_label, num_alleles));
        }

        return locusMap;
    }

    template <typename InfectionEvent, typename LocusImpl>
    std::vector<InfectionEvent*> parseInfectionsFromJSON(
            const json &input,
            std::map<std::string, LocusImpl*> loci,
            const char infectionsKey[] = "nodes",
            const char obsGenotypesKey[] = "observed_genotype",
            const char idKey[] = "id",
            const char observationDateKey[] = "observation_time",
            const char genotypeKey[] = "genotype",
            const char locusKey[] = "locus") {

        std::vector<InfectionEvent*> infections{};
        for (const auto& inf : input.at(infectionsKey)) {
            infections.push_back(new InfectionEvent(inf.at(idKey), inf.at(observationDateKey)));

            for (const auto& genetics : inf.at(obsGenotypesKey)) {
                auto locusLabel = genetics.at(locusKey);
                auto locusItr = loci.find(locusLabel);
                if(locusItr == loci.end()) {
                    std::cerr << "Locus " << locusLabel << " for node " << infections.back()->id() << " does not exist.\n";
                    exit(1);
                }
                if (!missingGenotype(genetics.at(genotypeKey))) {
                    std::string obs_genetics = genetics.at(genotypeKey);
                    infections.back()->addGenetics(locusItr->second, obs_genetics, obs_genetics);
                } else {
                    std::string latent_genetics;
                    latent_genetics.resize(locusItr->second->totalAlleles(), '0');
                    latent_genetics[0] = '1';
                    infections.back()->addLatentGenetics(locusItr->second, latent_genetics);
                }
            }
        }

        return infections;
    }

    template <typename InfectionEvent>
    std::map<InfectionEvent*, std::vector<InfectionEvent*>> parseDisallowedParentsFromJSON(
            const json &input,
            std::vector<InfectionEvent*> infections,
            const char infectionsKey[] = "nodes",
            const char idKey[] = "id",
            const char disallowedParentsKey[] = "disallowed_parents"
            ) {

        std::map<InfectionEvent*, std::vector<InfectionEvent*>> disallowedParents{};

        for(const auto& targetInfection : infections) {
            for (const auto& inf : input.at(infectionsKey)) {
                auto infectionId = inf.at(idKey);
                if (infectionId == targetInfection->id()) {
                    for (const auto& parentId : inf.at(disallowedParentsKey)) {
                        auto parentInf = std::find_if(infections.begin(), infections.end(), [&parentId](const InfectionEvent* candidateInf) {return candidateInf->id() == parentId;});
                        disallowedParents[targetInfection].push_back(*parentInf);
                    }
                    break;
                }
            }
        }

        return disallowedParents;

    }
}






#endif //TRANSMISSION_NETWORKS_APP_PARSE_JSON_H

