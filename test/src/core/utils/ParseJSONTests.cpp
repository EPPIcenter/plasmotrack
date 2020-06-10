//
// Created by Maxwell Murphy on 5/18/20.
//

#include <fstream>

#include "gtest/gtest.h"

#include "core/utils/io/parse_json.h"

#include "core/containers/Locus.h"
#include "core/containers/Infection.h"
#include "core/datatypes/Alleles.h"


TEST(ParseJSONTests, TestParseLoci) {
    using InfectionEvent = Infection<AllelesBitSet<32>>;
    std::ifstream testFile("/Users/maxwellmurphy/Workspace/transmission_nets/test/resources/JSON/nodes.json");
    if(!testFile) {
        std::cout << "Cannot open file." << std::endl;
        exit(1);
    }

    auto j = loadJSON(testFile);
    auto locusMap = parseLociFromJSON<Locus>(j);
    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, locusMap);

    ASSERT_TRUE(missingGenotype("0000000"));
    ASSERT_TRUE(missingGenotype(""));
    ASSERT_FALSE(missingGenotype("00001"));
}