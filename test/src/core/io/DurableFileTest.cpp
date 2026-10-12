//
// Writing a file durably: the whole new file replaces the old one under its name, nothing else is left behind, and a
// failure ends the process with OUTPUT_FAILURE naming the path.
//

#include "gtest/gtest.h"

#include "core/io/DurableFile.h"
#include "core/io/OutputFailure.h"

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

namespace fs = std::filesystem;
using transmission_nets::core::io::writeDurably;
using transmission_nets::core::io::writeFile;

namespace {
    std::string readFile(const fs::path& p) {
        std::ifstream in(p, std::ios::binary);
        return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    }

    struct DurableFileTest : ::testing::Test {
        fs::path dir = fs::temp_directory_path() / ("plasmotrack-durable-" + std::to_string(::getpid()) + "-" + ::testing::UnitTest::GetInstance()->current_test_info()->name());

        void SetUp() override { fs::remove_all(dir); }
        void TearDown() override {
            if (fs::exists(dir)) {
                fs::permissions(dir, fs::perms::owner_all);
            }
            if (!HasFailure()) {
                fs::remove_all(dir);
            }
        }
    };

    using DurableFileDeathTest = DurableFileTest;
}// namespace

TEST_F(DurableFileTest, WritesTheFileInItsDirectoryAndLeavesNothingElse) {
    const auto file = dir / "sub" / "checkpoint.cbor";
    writeDurably(file, [](const fs::path& temporary) { writeFile(temporary, "first"); });
    writeDurably(file, [](const fs::path& temporary) { writeFile(temporary, "second"); });

    EXPECT_EQ(readFile(file), "second");
    EXPECT_EQ(std::distance(fs::directory_iterator(dir / "sub"), fs::directory_iterator()), 1);
}

// The temporary file is complete and synced before it takes the name, so a failed rename leaves the old file whole.
TEST_F(DurableFileDeathTest, AFileThatCannotTakeItsNameExitsWithOutputFailureNamingThePath) {
    const auto file = dir / "checkpoint.cbor";
    fs::create_directories(file / "occupied");
    EXPECT_EXIT(writeDurably(file, [](const fs::path& temporary) { writeFile(temporary, "bytes"); }),
                ::testing::ExitedWithCode(transmission_nets::core::io::OUTPUT_FAILURE), "^ERROR: cannot rename into place output file " + file.string() + ": ");
}

TEST_F(DurableFileDeathTest, AnUnwritableDirectoryExitsWithOutputFailureNamingThePath) {
    fs::create_directories(dir);
    fs::permissions(dir, fs::perms::owner_read | fs::perms::owner_exec);
    EXPECT_EXIT(writeDurably(dir / "manifest.json", [](const fs::path& temporary) { writeFile(temporary, "{}"); }),
                ::testing::ExitedWithCode(transmission_nets::core::io::OUTPUT_FAILURE), "^ERROR: cannot open output file " + (dir / ".manifest.json.tmp").string() + ": ");
}
