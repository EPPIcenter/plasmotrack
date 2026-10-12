//
// Runs the CLI binary for the end-to-end tests and reads back what a Run left: its console output, exit code, Run events
// and closed Windows.
//

#ifndef TRANSMISSION_NETWORKS_TEST_CLIRUNNING_H
#define TRANSMISSION_NETWORKS_TEST_CLIRUNNING_H

#include "StoreReading.h"

#include "gtest/gtest.h"

#include <nlohmann/json.hpp>

#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <thread>
#include <vector>

extern char** environ;

namespace transmission_nets::test {
    inline fs::path makeWork(const std::string& name) {
        const fs::path work = fs::temp_directory_path() / ("plasmotrack-" + name + "-" + std::to_string(::getpid()));
        fs::remove_all(work);
        fs::create_directories(work / "out");
        return work;
    }

    // The CLI's arguments for the example input and the golden IDPs, then `out` as the output directory unless it is
    // empty, then `options`.
    inline std::vector<std::string> cliArguments(const fs::path& out, const std::vector<std::string>& options, const fs::path& input = EXAMPLE_INPUT_PATH) {
        const fs::path golden = GOLDEN_DATA_DIR;
        std::vector<std::string> args{GOLDEN_CLI_PATH, input.string()};
        if (!out.empty()) {
            args.push_back(out.string());
        }
        args.insert(args.end(), {"--symptomatic-idp", (golden / "symptomatic_idp.txt").string(), "--asymptomatic-idp", (golden / "asymptomatic_idp.txt").string()});
        args.insert(args.end(), options.begin(), options.end());
        return args;
    }

    // Starts the CLI with `args`, its standard output in <work>/cli.txt and its standard error in <work>/stderr.txt,
    // and returns its pid.
    inline pid_t spawnCli(const fs::path& work, std::vector<std::string> args) {
        std::vector<char*> argv;
        for (auto& a : args) {
            argv.push_back(a.data());
        }
        argv.push_back(nullptr);
        posix_spawn_file_actions_t actions;
        posix_spawn_file_actions_init(&actions);
        const auto outLog = (work / "cli.txt").string();
        const auto errLog = (work / "stderr.txt").string();
        posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, outLog.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, errLog.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        pid_t pid = -1;
        const int err = posix_spawn(&pid, GOLDEN_CLI_PATH, &actions, nullptr, argv.data(), environ);
        posix_spawn_file_actions_destroy(&actions);
        EXPECT_EQ(err, 0) << "cannot start the CLI";
        return pid;
    }

    // Starts the CLI writing a Run to <work>/out.
    inline pid_t startCli(const fs::path& work, const std::vector<std::string>& options) {
        return spawnCli(work, cliArguments(work / "out", options));
    }

    // The exit code, or 128 + the signal number for a CLI a signal ended, as a shell reports it.
    inline int waitFor(const pid_t pid) {
        int status = 0;
        ::waitpid(pid, &status, 0);
        return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    }

    inline int runCli(const fs::path& work, const std::vector<std::string>& options) {
        return waitFor(startCli(work, options));
    }

    inline std::string readText(const fs::path& file) {
        std::ifstream in(file);
        return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    }

    // What the last CLI started in `work` wrote to standard error.
    inline std::string errorOutput(const fs::path& work) {
        return readText(work / "stderr.txt");
    }

    // What the last CLI started in `work` printed: its standard output, then its standard error.
    inline std::string consoleOutput(const fs::path& work) {
        return readText(work / "cli.txt") + errorOutput(work);
    }

    // Waits until `condition` holds. False if the CLI exits or a minute passes first.
    template<typename Condition>
    bool waitUntil(const pid_t pid, Condition&& condition) {
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::minutes(1);
        while (!condition()) {
            if (std::chrono::steady_clock::now() > deadline or ::waitpid(pid, nullptr, WNOHANG) != 0) {
                return false;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        return true;
    }

    // Waits until the store has `count` closed Windows. False if the CLI exits or a minute passes first.
    inline bool waitForClosedWindows(const pid_t pid, const fs::path& store, const std::size_t count) {
        return waitUntil(pid, [&] { return closedWindows(store).size() >= count; });
    }

    // Waits until the last CLI started in `work` has printed `text`. False if it exits or a minute passes first.
    inline bool waitForConsoleOutput(const pid_t pid, const fs::path& work, const std::string& text) {
        return waitUntil(pid, [&] { return consoleOutput(work).find(text) != std::string::npos; });
    }

    inline nlohmann::json readJson(const fs::path& p) {
        std::ifstream in(p);
        return nlohmann::json::parse(in);
    }

    // Every closed Window's marker and Trace rows agree: each Window covers at least one Sampling iteration, each Trace
    // row falls in its Window's iterations and carries its Window number, and every table file of a closed Window reads
    // back.
    inline void expectEveryClosedWindowReadsBack(const fs::path& store) {
        const auto windows = closedWindows(store);
        ASSERT_FALSE(windows.empty());
        std::map<std::int64_t, nlohmann::json> markers;
        for (const auto& w : windows) {
            const auto marker            = readJson(store / "_closed" / w);
            markers[marker.at("window")] = marker;
            EXPECT_LE(marker.at("first_iteration").get<std::int64_t>(), marker.at("last_iteration").get<std::int64_t>()) << "Window " << w << " covers no Sampling iteration";
        }
        const auto tables = readStore(store);
        for (const auto& table : {"run_scalars", "infections", "ordering", "allele_frequencies", "latent_genotypes", "latent_parent_genotypes", "parent_sets", "chains"}) {
            const auto& rows = tables.at(table);
            ASSERT_GT(rows.rows(), 0u) << table;
            for (std::size_t r = 0; r < rows.rows(); ++r) {
                const auto window    = rows.column<std::int64_t>("window")[r];
                const auto iteration = rows.column<std::int64_t>("iteration")[r];
                ASSERT_TRUE(markers.contains(window)) << table << " has a row in Window " << window << ", which has no marker";
                EXPECT_GE(iteration, markers[window].at("first_iteration").get<std::int64_t>()) << table;
                EXPECT_LE(iteration, markers[window].at("last_iteration").get<std::int64_t>()) << table;
            }
        }
    }

    struct Event {
        std::int64_t iteration;
        std::string type;
        nlohmann::json payload;
    };

    inline std::vector<Event> events(const fs::path& store) {
        const auto t = readStoreTable(store, "events");
        std::vector<Event> out;
        for (std::size_t r = 0; r < t.rows(); ++r) {
            out.push_back({t.column<std::int64_t>("iteration")[r], t.column<std::string>("type")[r], nlohmann::json::parse(t.column<std::string>("payload")[r])});
        }
        return out;
    }

    inline std::vector<std::string> types(const std::vector<Event>& events) {
        std::vector<std::string> out;
        for (const auto& e : events) {
            out.push_back(e.type);
        }
        return out;
    }
}// namespace transmission_nets::test

#endif//TRANSMISSION_NETWORKS_TEST_CLIRUNNING_H
