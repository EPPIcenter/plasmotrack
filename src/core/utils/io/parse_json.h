//
// Created by Maxwell Murphy on 5/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARSE_JSON_H
#define TRANSMISSION_NETWORKS_APP_PARSE_JSON_H

#include <nlohmann/json.hpp>
#include <iostream>

#include "core/containers/Locus.h"
#include "core/containers/Infection.h"

using nlohmann::json;

inline json loadJSON(std::istream &input) {
    auto j = json::parse(input);
    return j;
};

inline bool missingGenotype(const std::string& alleleStr) {
    bool isEmpty = alleleStr.size() == 0;
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
        locusMap.emplace(locus_label, new Locus(locus_label, num_alleles));
    }

    return locusMap;
};

template <typename InfectionEvent, typename LocusImpl>
std::vector<InfectionEvent*> parseInfectionsFromJSON(
        const json &input,
        std::map<std::string, LocusImpl*> loci,
        const char infectionsKey[] = "nodes",
        const char obsGenotypesKey[] = "observed_genotype",
        const char idKey[] = "id",
        const char genotypeKey[] = "genotype",
        const char locusKey[] = "locus") {

    std::vector<InfectionEvent*> infections{};
    for (const auto& inf : input.at(infectionsKey)) {
        infections.push_back(new InfectionEvent(inf.at(idKey)));
        for (const auto& genetics : inf.at(obsGenotypesKey)) {
            auto locusLabel = genetics.at(locusKey);
            auto locusItr = loci.find(locusLabel);
            if(locusItr == loci.end()) {
                std::cerr << "Locus " << locusLabel << " for node " << infections.back()->id() << " does not exist.\n";
                exit(1);
            }
            if (!missingGenotype(genetics.at(genotypeKey))) {
                std::string obs_genetics = genetics.at(genotypeKey);
//                std::string latent_genetics = "";
//                latent_genetics.resize(obs_genetics.length(), '1');
                infections.back()->addGenetics(locusItr->second, obs_genetics, obs_genetics);
//                infections.back()->addGenetics(locusItr->second, obs_genetics, latent_genetics);
            }

        }
    }

    return infections;
}




#endif //TRANSMISSION_NETWORKS_APP_PARSE_JSON_H

