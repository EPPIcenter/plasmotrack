//
// The console progress adapter: the lines it writes for what the Run driver reports.
//

#include "gtest/gtest.h"

#include "impl/model/Model/Progress.h"

#include <sstream>
#include <string>

using namespace transmission_nets::impl::Model;

namespace {
    // Positions in a Run of 1000 burn-in and 1000 sampling iterations.
    RunPosition inBurnIn(const std::int64_t iteration, const std::int64_t lastClosedWindow) {
        return {.iteration = iteration, .phase = Phase::BurnIn, .lastClosedWindow = lastClosedWindow, .burnInRemaining = 1000 - iteration, .samplingRemaining = 1000};
    }
    RunPosition inSampling(const std::int64_t iteration, const std::int64_t lastClosedWindow) {
        return {.iteration = iteration, .phase = Phase::Sampling, .lastClosedWindow = lastClosedWindow, .burnInRemaining = 0, .samplingRemaining = 2000 - iteration};
    }
}// namespace

// Phase, iteration out of the total, the Window, the ETA at the Window's rate, the Window's and the Run's rates, the
// Hot chain's log prior and log likelihood, and swap acceptance per neighbouring pair of temperature slots, "-" for a
// pair with no swap attempted.
TEST(ConsoleProgressTest, AClosedWindowIsOneLine) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.windowClosed(inBurnIn(150, 2), {.window = 2, .phase = Phase::BurnIn, .samplingIterations = 50, .seconds = 4, .hotPrior = -1234.5678, .hotLikelihood = -5678.9, .swaps = {{.accepted = 21, .attempted = 50}, {}}});
    EXPECT_EQ(out.str(), "burn-in 150/2000  Window 2  ETA 2m28s  12.50 iterations/s (Run 12.50)  prior/likelihood -1234.57 / -5678.90  swap acceptance 0.42 -\n");
}

// The Run's rate is over every Window closed so far. One chain has no swap acceptance to show.
TEST(ConsoleProgressTest, TheRunsRateIsOverEveryWindowClosed) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.windowClosed(inSampling(1150, 11), {.window = 11, .phase = Phase::Sampling, .samplingIterations = 50, .seconds = 4, .hotPrior = -1, .hotLikelihood = -2});
    out.str("");
    console.windowClosed(inSampling(1250, 12), {.window = 12, .phase = Phase::Sampling, .samplingIterations = 100, .seconds = 4, .hotPrior = -1, .hotLikelihood = -2});
    EXPECT_EQ(out.str(), "sampling 1250/2000  Window 12  ETA 30s  25.00 iterations/s (Run 18.75)  prior/likelihood -1.00 / -2.00\n");
}

// The ETA shows hours and minutes once it is an hour or more; with no time spent, there is no rate to give.
TEST(ConsoleProgressTest, TheEtaIsInHoursAndMinutesOrUnknown) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.windowClosed(inBurnIn(100, 0), {.window = 0, .phase = Phase::BurnIn, .samplingIterations = 100, .seconds = 400, .hotPrior = -1, .hotLikelihood = -2});
    EXPECT_EQ(out.str(), "burn-in 100/2000  Window 0  ETA 2h06m  0.25 iterations/s (Run 0.25)  prior/likelihood -1.00 / -2.00\n");

    std::ostringstream still;
    ConsoleProgress stalled(still);
    stalled.windowClosed(inBurnIn(1, 0), {.window = 0, .phase = Phase::BurnIn, .samplingIterations = 1, .seconds = 0, .hotPrior = -1, .hotLikelihood = -2});
    EXPECT_EQ(still.str(), "burn-in 1/2000  Window 0  ETA ?  - iterations/s (Run -)  prior/likelihood -1.00 / -2.00\n");
}

// A heartbeat gives where the Run stands, the ETA at this process's rate so far and the Hot chain's log densities.
TEST(ConsoleProgressTest, AHeartbeatIsOneShortLine) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.heartbeat(inBurnIn(160, 2), {.hotPrior = -1, .hotLikelihood = -2, .samplingIterations = 60, .seconds = 120});
    EXPECT_EQ(out.str(), "burn-in 160/2000  heartbeat  ETA 1h01m  prior/likelihood -1.00 / -2.00\n");
}

// Ladder adaptations print no line of their own: the next Window line ends with the latest ladder and how many
// adaptations changed it since the last Window line. One whose temperatures came out NaN left the ladder as it was.
TEST(ConsoleProgressTest, LadderAdaptationsAreSummedUpOnTheNextWindowLine) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.ladderAdapted(inBurnIn(16, -1), {.applied = true, .temperatures = {1, 0.7877, 0.6667}, .barriers = {0.25, 0.75}});
    console.ladderAdapted(inBurnIn(21, -1), {.applied = false, .temperatures = {1, 0.7877, 0.6667}, .barriers = {0.25, 0.75}});
    console.ladderAdapted(inBurnIn(26, -1), {.applied = true, .temperatures = {1, 0.90821, 0.66667}, .barriers = {0.25, 0.75}});
    EXPECT_EQ(out.str(), "");
    const WindowSummary window{.window = 0, .phase = Phase::BurnIn, .samplingIterations = 100, .seconds = 4, .hotPrior = -1, .hotLikelihood = -2, .swaps = {{.accepted = 1, .attempted = 5}, {.accepted = 5, .attempted = 5}}};
    console.windowClosed(inBurnIn(100, 0), window);
    EXPECT_EQ(out.str(), "burn-in 100/2000  Window 0  ETA 1m16s  25.00 iterations/s (Run 25.00)  prior/likelihood -1.00 / -2.00  swap acceptance 0.20 1.00  ladder 1, 0.908, 0.667 (adapted 2×)\n");

    out.str("");
    console.ladderAdapted(inBurnIn(131, 0), {.applied = false, .temperatures = {1, 0.90821, 0.66667}, .barriers = {0.25, 0.75}});
    console.windowClosed(inBurnIn(200, 1), {.window = 1, .phase = Phase::BurnIn, .samplingIterations = 100, .seconds = 4, .hotPrior = -1, .hotLikelihood = -2, .swaps = window.swaps});
    EXPECT_EQ(out.str(), "burn-in 200/2000  Window 1  ETA 1m12s  25.00 iterations/s (Run 25.00)  prior/likelihood -1.00 / -2.00  swap acceptance 0.20 1.00\n");
}

TEST(ConsoleProgressTest, TheRunsStartAndEndAreOneLineEach) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.started({.iteration = 0, .phase = Phase::BurnIn, .lastClosedWindow = -1, .burnInRemaining = 1000, .samplingRemaining = 1000}, false);
    console.started(inSampling(1250, 12), true);
    console.stopped(inSampling(1300, 13), "SIGINT");
    console.ended(inSampling(2000, 20));
    EXPECT_EQ(out.str(), "Starting at iteration 0/2000.\n"
                         "Resuming at iteration 1250/2000 (Window 12 committed).\n"
                         "Stopped at iteration 1300 (Window 13 committed).\n"
                         "Completed at iteration 2000 (Window 20 committed).\n");
}

// A notice before the Run starts is its own line, as given.
TEST(ConsoleProgressTest, ANoticeIsOneLine) {
    std::ostringstream out;
    ConsoleProgress console(out);
    console.notice("Seed Used: 7");
    EXPECT_EQ(out.str(), "Seed Used: 7\n");
}
