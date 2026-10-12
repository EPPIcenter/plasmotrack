//
// Golden run: the CLI, run on the example input with a fixed seed, must write the same store as the golden store under
// test/data/golden/<variant>/: the same files, the same Window markers, the same decoded tables and the same Run
// manifest. What differs from one Run to the next by design is left out: the manifest's Run id, start time and commit,
// the paths in its CLI arguments, and every wall time and duration in the events table. Comparing decoded tables also
// leaves out Parquet file metadata such as `created_by`.
//
// The golden stores were written by these commands (the example IDPs give prior -inf, so working IDPs are committed
// next to them), and their manifests then put through `normalised` below:
//
//   transmission_networks_model examples/example_input.json OUT \
//       --symptomatic-idp test/data/golden/symptomatic_idp.txt --asymptomatic-idp test/data/golden/asymptomatic_idp.txt \
//       --seed 42 -b 20 -s 50 -t 5 --window 25 [-n 4 -g 0.5 -c 1]
//
// Before they replaced it, the store of these commands (with the default --window) held the same values as the gzip
// capture of #39 for every quantity both hold.
//
// Never rewrite a golden store to make this pass after a change: a mismatch means the output changed.
//

#include "gtest/gtest.h"

#include "StoreReading.h"

#include <nlohmann/json.hpp>

#include <unistd.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using transmission_nets::core::io::Table;
using transmission_nets::test::readStore;
using transmission_nets::test::withoutWallTimes;

namespace {
    std::string readFile(const fs::path& p) {
        std::ifstream in(p, std::ios::binary);
        return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    }

    std::set<std::string> filesUnder(const fs::path& root) {
        std::set<std::string> out;
        for (const auto& e : fs::recursive_directory_iterator(root)) {
            if (e.is_regular_file()) out.insert(fs::relative(e.path(), root).generic_string());
        }
        return out;
    }

    // The Run manifest without what differs between Runs of one recipe; the paths in the CLI arguments become names.
    nlohmann::json normalised(nlohmann::json manifest, const std::map<std::string, std::string>& paths) {
        manifest["run_id"]     = "<run id>";
        manifest["started_at"] = "<wall time>";
        manifest["git_commit"] = "<commit>";
        for (auto& argument : manifest.at("cli_arguments")) {
            auto text = argument.get<std::string>();
            for (const auto& [path, name] : paths) {
                if (const auto at = text.find(path); at != std::string::npos) {
                    text.replace(at, path.size(), name);
                }
            }
            argument = text;
        }
        return manifest;
    }

    void expectSameTable(const std::string& name, const Table& actual, const Table& expected) {
        ASSERT_EQ(actual.names(), expected.names()) << name;
        ASSERT_EQ(actual.rows(), expected.rows()) << name;
        for (std::size_t c = 0; c < expected.columns.size(); ++c) {
            EXPECT_EQ(actual.columns[c].second, expected.columns[c].second) << name << "." << expected.columns[c].first;
        }
    }

    struct GoldenVariant {
        std::string name;
        std::string extraArgs;

        friend std::ostream& operator<<(std::ostream& os, const GoldenVariant& v) { return os << v.name; }
    };

    class GoldenRunTest : public ::testing::TestWithParam<GoldenVariant> {};
}// namespace

TEST_P(GoldenRunTest, StoreMatchesTheGoldenStore) {
    const auto& variant   = GetParam();
    const fs::path golden = GOLDEN_DATA_DIR;
    const fs::path expect = golden / variant.name;
    ASSERT_TRUE(fs::exists(expect / "manifest.json")) << "missing golden store " << expect;

    const fs::path work = fs::temp_directory_path() / ("plasmotrack-golden-" + variant.name + "-" + std::to_string(::getpid()));
    fs::remove_all(work);
    const fs::path out = work / "out";
    fs::create_directories(out);

    std::ostringstream cmd;
    cmd << '"' << GOLDEN_CLI_PATH << "\" \"" << EXAMPLE_INPUT_PATH << "\" \"" << out.string() << '"'
        << " --symptomatic-idp \"" << (golden / "symptomatic_idp.txt").string() << '"'
        << " --asymptomatic-idp \"" << (golden / "asymptomatic_idp.txt").string() << '"'
        << " --seed 42 -b 20 -s 50 -t 5 --window 25 " << variant.extraArgs
        << " > \"" << (work / "cli.txt").string() << "\" 2>&1";
    ASSERT_EQ(std::system(cmd.str().c_str()), 0) << "CLI failed; console output in " << work / "cli.txt";

    // Checkpoints carry the Run id, so they are not part of the golden store; their contents are tested elsewhere.
    auto files = filesUnder(out);
    for (const auto* checkpoint : {"checkpoint.cbor", "checkpoint-burn-in.cbor"}) {
        EXPECT_EQ(files.erase(checkpoint), 1u) << checkpoint;
    }
    EXPECT_EQ(files, filesUnder(expect)) << "the set of files in the store changed";
    for (const auto& file : filesUnder(expect)) {
        if (file.starts_with("_closed/")) {
            EXPECT_EQ(nlohmann::json::parse(readFile(out / file)), nlohmann::json::parse(readFile(expect / file))) << file;
        }
    }

    const std::map<std::string, std::string> paths{{GOLDEN_CLI_PATH, "<cli>"}, {EXAMPLE_INPUT_PATH, "<input>"}, {out.string(), "<out>"}, {golden.string(), "<golden>"}};
    EXPECT_EQ(normalised(nlohmann::json::parse(readFile(out / "manifest.json")), paths), nlohmann::json::parse(readFile(expect / "manifest.json")));

    const auto actual   = readStore(out);
    const auto expected = readStore(expect);
    ASSERT_EQ(actual.size(), expected.size());
    for (const auto& [name, table] : expected) {
        ASSERT_TRUE(actual.contains(name)) << name;
        if (name == "events") {
            expectSameTable(name, withoutWallTimes(actual.at(name)), withoutWallTimes(table));
        } else {
            expectSameTable(name, actual.at(name), table);
        }
    }

    if (!HasFailure()) fs::remove_all(work);
}

INSTANTIATE_TEST_SUITE_P(
        CliOutput, GoldenRunTest,
        ::testing::Values(
                GoldenVariant{"single-chain", ""},
                GoldenVariant{"replica-exchange", "-n 4 -g 0.5 -c 1"}),
        [](const auto& info) {
            auto n = info.param.name;
            std::erase(n, '-');
            return n;
        });
