//
// The JSON adapter of the Study. See Study.h.
//

#include "Study.h"

#include "Refusal.h"

#include <fmt/core.h>

namespace transmission_nets::impl::Model {
    Study studyFromJson(const nlohmann::json& input, std::vector<double> symptomaticInfectionDurations, std::vector<double> asymptomaticInfectionDurations) {
        Study study{.symptomaticInfectionDurations = std::move(symptomaticInfectionDurations), .asymptomaticInfectionDurations = std::move(asymptomaticInfectionDurations)};
        try {
            for (const auto& locus : input.at("loci")) {
                study.loci.push_back({.label = locus.at("locus"), .alleles = locus.at("num_alleles"), .alleleFrequencies = locus.at("allele_freqs")});
            }
            for (const auto& node : input.at("nodes")) {
                Study::Infection infection{.id = node.at("id"), .observationTime = node.at("observation_time")};
                if (node.contains("symptomatic")) {
                    infection.symptomatic = node.at("symptomatic");
                }
                for (const auto& genotype : node.at("observed_genotype")) {
                    infection.observedGenotypes.push_back({.locus = genotype.at("locus"), .alleles = genotype.at("genotype")});
                }
                // R's jsonlite writes an empty list as {}.
                const auto& allowedParents = node.at("allowed_parents");
                if (!(allowedParents.is_object() && allowedParents.empty())) {
                    infection.permittedParents = allowedParents.get<std::vector<std::string>>();
                }
                study.infections.push_back(std::move(infection));
            }
        } catch (const nlohmann::json::exception& e) {
            throw Refused(refuse(RefusalKind::InvalidInput, fmt::format("the input does not describe a study: {}", e.what())));
        }
        return study;
    }
}// namespace transmission_nets::impl::Model
