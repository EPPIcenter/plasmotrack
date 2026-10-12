//
// How the CLI fails: an output it cannot write ends the Run with exit code 3 and a stderr message naming the path,
// instead of exiting 0; a refused Run exits 1 with the refusal on stderr, and writes nothing.
//

#include "gtest/gtest.h"

#include "CliRunning.h"

#include <filesystem>
#include <string>
#include <vector>

using namespace transmission_nets::test;

namespace {
    // A short Run with `options` added.
    std::vector<std::string> shortRun(const std::vector<std::string>& options = {"-t", "5"}) {
        std::vector<std::string> out{"--seed", "42", "-b", "2", "-s", "5"};
        out.insert(out.end(), options.begin(), options.end());
        return out;
    }
}// namespace

TEST(CliFailure, UnwritableOutputDirExitsNonZeroNamingThePath) {
    const auto work = makeWork("cli-failure");
    const auto out  = work / "out";
    fs::permissions(out, fs::perms::owner_read | fs::perms::owner_exec);

    const int exitCode = runCli(work, shortRun());

    fs::permissions(out, fs::perms::owner_all);
    EXPECT_EQ(exitCode, 3) << "expected the output-failure exit code, not 0 or a crash";
    EXPECT_NE(errorOutput(work).find(out.string()), std::string::npos) << "stderr: " << errorOutput(work);

    if (!HasFailure()) fs::remove_all(work);
}

// Which options a Run refuses is tested in RunRecipeTest; this checks that a refusal reaches stderr, exits 1 and
// writes nothing.
TEST(CliFailure, ARefusalReachesStderrExitsOneAndWritesNothing) {
    const auto work = makeWork("cli-refusal");

    EXPECT_EQ(runCli(work, shortRun({"-t", "5", "--window", "12"})), 1);
    EXPECT_NE(errorOutput(work).find("ERROR: --window must be a positive multiple of -t"), std::string::npos) << "stderr: " << errorOutput(work);
    EXPECT_TRUE(fs::is_empty(work / "out")) << "a refused Run writes no output";

    if (!HasFailure()) fs::remove_all(work);
}
