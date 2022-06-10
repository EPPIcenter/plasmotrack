//
// Created by Maxwell Murphy on 5/18/20.
//

#include <fstream>

#include "gtest/gtest.h"

#include "core/io/parse_json.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Alleles.h"

using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::io;

TEST(ParseJSONTests, TestMissingGenotype) {
    ASSERT_TRUE(missingGenotype("0000000"));
    ASSERT_TRUE(missingGenotype(""));
    ASSERT_FALSE(missingGenotype("00001"));
}