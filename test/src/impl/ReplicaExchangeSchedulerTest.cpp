//
// Replica exchange builds each chain's scheduler from (State, Model, rng) alone. Every scheduler the Model offers must
// be constructible that way; this test instantiates ReplicaExchange with each, so a scheduler that drifts from that
// shape stops the test build rather than going unnoticed until someone switches the CLI to it.
//

#include "gtest/gtest.h"

#include "core/samplers/meta/ReplicaExchange.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/SampleScheduler.h"
#include "impl/model/Model/SequentialScheduler.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <memory>
#include <optional>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::samplers::ReplicaExchange;

namespace {
    template<template<typename...> typename Scheduler>
    void runTwoChains() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        const auto input = nlohmann::json::parse(in);
        const std::vector<double> idp(60, 1.0 / 60);
        auto rng         = std::make_shared<boost::random::mt19937>(20261011);

        ReplicaExchange<State, Model, Scheduler> repex(2, 0.5, rng, false, 1, input, idp, idp);
        repex.burnin();
        repex.sample();

        EXPECT_TRUE(std::isfinite(repex.hotValue()));
    }
}// namespace

TEST(ReplicaExchangeSchedulerTest, RunsWithTheSequentialSampleScheduler) {
    runTwoChains<SequentialSampleScheduler>();
}

TEST(ReplicaExchangeSchedulerTest, RunsWithTheSampleScheduler) {
    runTwoChains<SampleScheduler>();
}
