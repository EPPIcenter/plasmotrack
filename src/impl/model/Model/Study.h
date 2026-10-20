//
// The Study: the input as the model reads it. Loci with their Allele counts and frequencies; Infections with their
// observation time, whether symptomatic, their Observed genotypes and the parents the study permits them; and the
// infection-duration distributions for symptomatic and asymptomatic Infections. Everything is in input order. It is
// read once from the input files, by studyFromJson, and a State is built from it; a test builds one as a literal. The
// rules the model applies to the data (starting genotypes, the frequency floor, the null model) are the State's.
//

#ifndef TRANSMISSION_NETWORKS_APP_STUDY_H
#define TRANSMISSION_NETWORKS_APP_STUDY_H

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace transmission_nets::impl::Model {
    struct Study {
        struct Locus {
            std::string label;
            unsigned int alleles = 0;
            std::vector<double> alleleFrequencies;// one per Allele; the State floors them
        };

        // An Observed genotype at one Locus: one '0' or '1' per Allele. All zeros means nothing was observed there, as
        // does a Locus with no entry.
        struct Genotype {
            std::string locus;
            std::string alleles;
        };

        struct Infection {
            std::string id;
            double observationTime = 0;
            bool symptomatic       = true;
            std::vector<Genotype> observedGenotypes;
            std::vector<std::string> permittedParents;// ids of the Infections the study allows as this one's parents
        };

        std::vector<Locus> loci;
        std::vector<Infection> infections;
        // Probability of each infection duration, by day since infection, for symptomatic and for asymptomatic
        // Infections.
        std::vector<double> symptomaticInfectionDurations;
        std::vector<double> asymptomaticInfectionDurations;
    };

    // The JSON adapter: the nodes file (`loci`: `locus`, `num_alleles`, `allele_freqs`; `nodes`: `id`,
    // `observation_time`, `symptomatic` (default true), `observed_genotype` of `locus` and `genotype`,
    // `allowed_parents`) and the two IDP files as a Study. Throws Refused for a document that lacks a field or has
    // one of the wrong type.
    Study studyFromJson(const nlohmann::json& input, std::vector<double> symptomaticInfectionDurations, std::vector<double> asymptomaticInfectionDurations);
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_STUDY_H
