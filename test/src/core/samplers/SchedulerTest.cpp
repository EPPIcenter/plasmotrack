//
// Samplers adapt only inside their adaptation window and only when the step allows it, so that adaptation can stop
// when burn-in ends.
//

#include "gtest/gtest.h"

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/scheduler/RandomizedScheduler.h"
#include "core/samplers/scheduler/Scheduler.h"

#include <fmt/format.h>

#include <memory>
#include <stdexcept>

using namespace transmission_nets::core::samplers;

namespace {
    struct CountingSampler : AbstractSampler {
        int* updates;
        int* adaptations;

        CountingSampler(int* updates, int* adaptations) : updates(updates), adaptations(adaptations) {}
        void update() noexcept override { ++*updates; }
        void adapt() noexcept override { ++*adaptations; }
    };
}// namespace

TEST(SchedulerTest, AdaptsInsideTheWindowOnlyWhenAllowed) {
    int updates     = 0;
    int adaptations = 0;
    Scheduler scheduler;
    scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .adaptationStart = 2, .adaptationEnd = 6});

    for (int s = 0; s < 4; ++s) {
        scheduler.step(true);// steps 0-3: 2 and 3 are in the window
    }
    EXPECT_EQ(adaptations, 2);
    for (int s = 0; s < 4; ++s) {
        scheduler.step(false);// steps 4-7: 4 and 5 are in the window, but adaptation is off
    }
    EXPECT_EQ(adaptations, 2);
    EXPECT_EQ(updates, 8);
}

TEST(RandomizedSchedulerTest, MakesItsUpdatesPerIterationEachIteration) {
    int updates     = 0;
    int adaptations = 0;
    RandomizedScheduler<> scheduler(std::make_shared<boost::random::mt19937>(1), 3);
    scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = "CountingSampler/a"});

    scheduler.step();
    scheduler.step();
    EXPECT_EQ(updates, 6);
}

TEST(RandomizedSchedulerTest, WithoutACountMakesOneUpdatePerRegisteredSamplerEachIteration) {
    int updates     = 0;
    int adaptations = 0;
    RandomizedScheduler<> scheduler(std::make_shared<boost::random::mt19937>(1));
    for (int s = 0; s < 3; ++s) {
        scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = fmt::format("CountingSampler/{}", s)});
    }

    scheduler.step();
    scheduler.step();
    EXPECT_EQ(updates, 6);
}

TEST(SchedulerTest, RefusesASecondSamplerWithTheSameId) {
    int updates     = 0;
    int adaptations = 0;
    Scheduler scheduler;
    scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = "CountingSampler/a"});
    EXPECT_THROW(scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = "CountingSampler/a"}), std::invalid_argument);
}

TEST(RandomizedSchedulerTest, RefusesASecondSamplerWithTheSameId) {
    int updates     = 0;
    int adaptations = 0;
    RandomizedScheduler<> scheduler(std::make_shared<boost::random::mt19937>(1));
    scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = "CountingSampler/a"});
    EXPECT_THROW(scheduler.registerSampler({.sampler = std::make_unique<CountingSampler>(&updates, &adaptations), .id = "CountingSampler/a"}), std::invalid_argument);
}
