//
// Studies for tests: the example input as a Study, and the pieces a Study literal is built from.
//

#ifndef TRANSMISSION_NETWORKS_TEST_STUDIES_H
#define TRANSMISSION_NETWORKS_TEST_STUDIES_H

#include "impl/model/Model/Study.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <string>
#include <vector>

namespace transmission_nets::test::studies {
    using transmission_nets::impl::Model::Study;

    // A flat infection-duration distribution over `days` days.
    inline std::vector<double> flatDurations(const int days = 60) {
        return std::vector<double>(static_cast<size_t>(days), 1.0 / days);
    }

    // A Locus whose Alleles are equally frequent.
    inline Study::Locus uniformLocus(const std::string& label, const unsigned int alleles) {
        return {.label = label, .alleles = alleles, .alleleFrequencies = std::vector<double>(alleles, 1.0 / alleles)};
    }

    // The example input (Infections 1, 2 and 3; Loci AS1 and AS11) with flat infection-duration distributions.
    inline Study exampleStudy(const std::vector<double>& durations = flatDurations()) {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return impl::Model::studyFromJson(nlohmann::json::parse(in), durations, durations);
    }
}// namespace transmission_nets::test::studies

#endif//TRANSMISSION_NETWORKS_TEST_STUDIES_H
