//
// The Parquet store sink, driven through the StoreSink interface and read back with the Arrow C++ reader.
//

#include "gtest/gtest.h"

#include "core/io/ParquetStoreSink.h"

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/column_page.h>
#include <parquet/column_reader.h>

#include <unistd.h>

#include <algorithm>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;
using transmission_nets::core::io::ParquetStoreSink;
using transmission_nets::core::io::Table;
using transmission_nets::core::io::WindowMarker;

namespace {
    using I = std::vector<std::int64_t>;
    using D = std::vector<double>;
    using S = std::vector<std::string>;
    using B = std::vector<bool>;

    std::unique_ptr<parquet::arrow::FileReader> openParquet(const fs::path& file) {
        auto in = arrow::io::ReadableFile::Open(file.string()).ValueOrDie();
        return parquet::arrow::OpenFile(in, arrow::default_memory_pool()).ValueOrDie();
    }

    std::shared_ptr<arrow::Table> readParquet(const fs::path& file) {
        std::shared_ptr<arrow::Table> out;
        EXPECT_TRUE(openParquet(file)->ReadTable(&out).ok());
        return out->CombineChunks().ValueOrDie();
    }

    // The reader rule: a Window exists only if its marker exists. Returns a table's rows over every closed Window.
    std::vector<std::shared_ptr<arrow::Table>> readClosed(const fs::path& root, const std::string& table) {
        std::vector<std::string> windows;
        if (fs::exists(root / "_closed")) {
            for (const auto& marker : fs::directory_iterator(root / "_closed")) {
                if (const auto name = marker.path().filename().string(); !name.starts_with(".")) {
                    windows.push_back(name);
                }
            }
        }
        std::ranges::sort(windows);
        std::vector<std::shared_ptr<arrow::Table>> out;
        for (const auto& w : windows) {
            if (const auto file = root / table / (w + ".parquet"); fs::exists(file)) {
                out.push_back(readParquet(file));
            }
        }
        return out;
    }

    template<typename ArrayT>
    std::shared_ptr<ArrayT> column(const arrow::Table& t, const std::string& name) {
        return std::static_pointer_cast<ArrayT>(t.GetColumnByName(name)->chunk(0));
    }

    std::vector<std::string> filesIn(const fs::path& dir) {
        std::vector<std::string> out;
        if (fs::exists(dir)) {
            for (const auto& f : fs::directory_iterator(dir)) {
                out.push_back(f.path().filename().string());
            }
        }
        std::ranges::sort(out);
        return out;
    }

    Table infections(const std::int64_t window, const std::string& phase) {
        return {.columns = {{"iteration", I{0, 10}},
                            {"phase", S{phase, phase}},
                            {"window", I{window, window}},
                            {"infection", S{"inf_a", "inf_b"}},
                            {"expected_false_positives", D{0.25, 1e-12}},
                            {"symptomatic", B{true, false}}}};
    }

    class ParquetStoreSinkTest : public ::testing::Test {
    protected:
        void SetUp() override {
            root = fs::temp_directory_path() / ("plasmotrack-parquet-store-" + std::to_string(::getpid()) + "-" +
                                                ::testing::UnitTest::GetInstance()->current_test_info()->name());
            fs::remove_all(root);
        }

        void TearDown() override {
            if (fs::exists(root)) {
                fs::permissions(root, fs::perms::owner_all, fs::perm_options::add);
            }
            fs::remove_all(root);
        }

        fs::path root;
    };

    // gtest runs death tests first, before reading Parquet back has started Arrow's thread pool and made fork()
    // unsafe. The kill test reads back, so it comes last.
    using ParquetStoreSinkDeathTest = ParquetStoreSinkTest;
}// namespace

TEST_F(ParquetStoreSinkTest, CommittedWindowReadsBackWithItsSchemaAndValues) {
    {
        ParquetStoreSink sink(root);
        sink.appendRows("infections", infections(0, "burn_in"));
        sink.commitWindow({.window = 0, .firstIteration = 0, .lastIteration = 10, .phase = "burn_in"});
    }

    ASSERT_EQ(filesIn(root / "infections"), std::vector<std::string>{"window-000000.parquet"});
    const auto read = readParquet(root / "infections" / "window-000000.parquet");
    const auto expected = arrow::schema({arrow::field("iteration", arrow::int64()),
                                         arrow::field("phase", arrow::utf8()),
                                         arrow::field("window", arrow::int64()),
                                         arrow::field("infection", arrow::utf8()),
                                         arrow::field("expected_false_positives", arrow::float64()),
                                         arrow::field("symptomatic", arrow::boolean())});
    EXPECT_TRUE(read->schema()->Equals(*expected)) << read->schema()->ToString();
    ASSERT_EQ(read->num_rows(), 2);
    EXPECT_EQ(column<arrow::Int64Array>(*read, "iteration")->Value(1), 10);
    EXPECT_EQ(column<arrow::StringArray>(*read, "phase")->GetString(0), "burn_in");
    EXPECT_EQ(column<arrow::StringArray>(*read, "infection")->GetString(1), "inf_b");
    EXPECT_EQ(column<arrow::DoubleArray>(*read, "expected_false_positives")->Value(1), 1e-12);
    EXPECT_TRUE(column<arrow::BooleanArray>(*read, "symptomatic")->Value(0));
    EXPECT_FALSE(column<arrow::BooleanArray>(*read, "symptomatic")->Value(1));
}

TEST_F(ParquetStoreSinkTest, HandedOverTablesStayHiddenUntilTheWindowCommitsThenTheMarkerRecordsIt) {
    ParquetStoreSink sink(root);
    sink.appendRows("infections", infections(3, "sampling"));

    EXPECT_TRUE(readClosed(root, "infections").empty());
    for (const auto& name : filesIn(root / "infections")) {
        EXPECT_TRUE(name.starts_with(".")) << name;
    }
    EXPECT_TRUE(filesIn(root / "_closed").empty());

    sink.commitWindow({.window = 3, .firstIteration = 21, .lastIteration = 30, .phase = "sampling"});

    EXPECT_EQ(filesIn(root / "infections"), std::vector<std::string>{"window-000003.parquet"});
    ASSERT_EQ(filesIn(root / "_closed"), std::vector<std::string>{"window-000003"});
    std::ifstream marker(root / "_closed" / "window-000003");
    EXPECT_EQ(nlohmann::json::parse(marker),
              (nlohmann::json{{"window", 3}, {"first_iteration", 21}, {"last_iteration", 30}, {"phase", "sampling"}}));
    EXPECT_EQ(readClosed(root, "infections").size(), 1);
}

// #56: a Window's rows come a few at a time while it is open. A table's row group ends at the first point the caller
// allows once its rows reach the sink's budget, and goes to the table's one in-flight file at once; a table that stays
// under budget writes one row group per Window. The Window still commits as one file per table, read back in the order
// handed over.
TEST_F(ParquetStoreSinkTest, RowsOfAWindowCommitAsOneFilePerTableInRowGroupsCutAtTheBudget) {
    const auto rows = [](const std::int64_t iteration, const std::string& infection) {
        return Table{.columns = {{"iteration", I{iteration, iteration}}, {"phase", S{"sampling", "sampling"}}, {"window", I{2, 2}}, {"infection", S{infection, infection}}}};
    };
    // In Arrow, two rows of `rows` take 2 * (8 + (4 + 8) + 8 + (4 + 5)) = 74 bytes: two calls pass a budget of 100.
    ParquetStoreSink sink(root, 100);
    sink.appendRows("infections", rows(11, "inf_a"));
    sink.appendRows("chains", rows(11, "inf_a"));
    sink.allowRowGroupEnd();
    sink.appendRows("infections", rows(12, "inf_b"));
    sink.allowRowGroupEnd();
    sink.appendRows("infections", rows(13, "inf_c"));
    sink.allowRowGroupEnd();
    sink.appendRows("infections", rows(14, "inf_d"));

    for (const auto* table : {"infections", "chains"}) {
        ASSERT_EQ(filesIn(root / table).size(), 1u) << table;
        EXPECT_TRUE(filesIn(root / table).front().starts_with(".")) << table;
    }
    EXPECT_TRUE(filesIn(root / "_closed").empty());

    sink.commitWindow({.window = 2, .firstIteration = 11, .lastIteration = 14, .phase = "sampling"});

    EXPECT_EQ(filesIn(root / "infections"), std::vector<std::string>{"window-000002.parquet"});
    EXPECT_EQ(filesIn(root / "chains"), std::vector<std::string>{"window-000002.parquet"});
    const auto groups = [&](const std::string& table) {
        const auto meta = openParquet(root / table / "window-000002.parquet")->parquet_reader()->metadata();
        std::vector<std::int64_t> out;
        for (int g = 0; g < meta->num_row_groups(); ++g) {
            out.push_back(meta->RowGroup(g)->num_rows());
        }
        return out;
    };
    EXPECT_EQ(groups("infections"), (std::vector<std::int64_t>{4, 4}));
    EXPECT_EQ(groups("chains"), (std::vector<std::int64_t>{2}));
    const auto read = readParquet(root / "infections" / "window-000002.parquet");
    ASSERT_EQ(read->num_rows(), 8);
    for (std::int64_t r = 0; r < 8; ++r) {
        EXPECT_EQ(column<arrow::Int64Array>(*read, "iteration")->Value(r), 11 + r / 2);
        EXPECT_EQ(column<arrow::StringArray>(*read, "infection")->GetString(r), std::string("inf_") + static_cast<char>('a' + r / 2));
    }
}

// Window numbers are zero-padded to six digits and grow past them rather than being cut.
TEST_F(ParquetStoreSinkTest, WindowNamesHaveAtLeastSixDigitsAndAreNeverTruncated) {
    ParquetStoreSink sink(root);
    sink.appendRows("infections", infections(7, "sampling"));
    sink.commitWindow({.window = 7, .firstIteration = 1, .lastIteration = 10, .phase = "sampling"});
    sink.appendRows("infections", infections(1234567, "sampling"));
    sink.commitWindow({.window = 1234567, .firstIteration = 11, .lastIteration = 20, .phase = "sampling"});

    EXPECT_EQ(filesIn(root / "infections"), (std::vector<std::string>{"window-000007.parquet", "window-1234567.parquet"}));
    EXPECT_EQ(filesIn(root / "_closed"), (std::vector<std::string>{"window-000007", "window-1234567"}));
}

// Resuming removes every Window after the checkpoint's, whatever the width of its number, and leaves files that aren't
// Windows alone.
TEST_F(ParquetStoreSinkTest, DiscardingWindowsReadsNumbersOfAnyWidthAndIgnoresOtherFiles) {
    {
        ParquetStoreSink sink(root);
        for (const std::int64_t window : {999999, 1000000, 1000001}) {
            sink.appendRows("infections", infections(window, "sampling"));
            sink.commitWindow({.window = window, .firstIteration = window, .lastIteration = window, .phase = "sampling"});
        }
    }
    for (const auto* stray : {"w", "window-", "window-12x.parquet", "notes.txt"}) {
        std::ofstream(root / "infections" / stray) << "stray";
        std::ofstream(root / "_closed" / stray) << "stray";
    }

    ParquetStoreSink(root).discardWindowsAfter(1000000);

    EXPECT_EQ(filesIn(root / "infections"), (std::vector<std::string>{"notes.txt", "w", "window-", "window-1000000.parquet", "window-12x.parquet", "window-999999.parquet"}));
    EXPECT_EQ(filesIn(root / "_closed"), (std::vector<std::string>{"notes.txt", "w", "window-", "window-1000000", "window-12x.parquet", "window-999999"}));
}

// A Run killed after Window 3's marker but before its checkpoint, with Window 4 in flight, resumes from the checkpoint
// taken after Window 1: Windows 2 and 3 and every in-flight file go, markers and tables alike, and nothing else does.
TEST_F(ParquetStoreSinkTest, DiscardingWindowsRemovesLaterWindowsAndEveryFileInFlight) {
    {
        ParquetStoreSink sink(root);
        sink.writeStructure("loci", infections(0, "burn_in"));
        for (const std::int64_t window : {0, 1, 2, 3}) {
            sink.appendRows("infections", infections(window, "sampling"));
            sink.commitWindow({.window = window, .firstIteration = window, .lastIteration = window, .phase = "sampling"});
        }
        sink.appendRows("infections", infections(4, "sampling"));
    }
    ASSERT_EQ(filesIn(root / "infections").size(), 5u) << "Window 4's table is not left in flight";

    ParquetStoreSink(root).discardWindowsAfter(1);

    EXPECT_EQ(filesIn(root / "infections"), (std::vector<std::string>{"window-000000.parquet", "window-000001.parquet"}));
    EXPECT_EQ(filesIn(root / "_closed"), (std::vector<std::string>{"window-000000", "window-000001"}));
    EXPECT_EQ(filesIn(root / "structure"), (std::vector<std::string>{"loci.parquet"}));
}

TEST_F(ParquetStoreSinkTest, ManifestAndStructureAreWrittenAtRunStart) {
    ParquetStoreSink sink(root);
    const nlohmann::json manifest{{"run_id", "r1"}, {"seed", 42}, {"parent_set_floor", 1e-9}};
    sink.writeManifest(manifest);
    sink.writeStructure("loci", {.columns = {{"locus", S{"L1", "L1", "L2"}}, {"allele", I{0, 1, 0}}}});

    std::ifstream in(root / "manifest.json");
    EXPECT_EQ(nlohmann::json::parse(in), manifest);
    EXPECT_EQ(filesIn(root / "structure"), std::vector<std::string>{"loci.parquet"});
    const auto loci = readParquet(root / "structure" / "loci.parquet");
    ASSERT_EQ(loci->num_rows(), 3);
    EXPECT_EQ(column<arrow::StringArray>(*loci, "locus")->GetString(2), "L2");
    EXPECT_EQ(column<arrow::Int64Array>(*loci, "allele")->Value(1), 1);
    EXPECT_EQ(filesIn(root), (std::vector<std::string>{"manifest.json", "structure"}));
}

TEST_F(ParquetStoreSinkTest, TableWithNoRowsInAWindowHasNoFile) {
    ParquetStoreSink sink(root);
    Table empty = infections(0, "burn_in");
    for (auto& [_, c] : empty.columns) {
        std::visit([](auto& v) { v.clear(); }, c);
    }
    sink.appendRows("infections", empty);
    sink.commitWindow({.window = 0, .firstIteration = 0, .lastIteration = 10, .phase = "burn_in"});

    EXPECT_FALSE(fs::exists(root / "infections"));
    EXPECT_EQ(filesIn(root / "_closed"), std::vector<std::string>{"window-000000"});
}

// #37: dictionary encoding off for `parent_sets.set_id` only, pages of 2^20 rows, zstd everywhere.
TEST_F(ParquetStoreSinkTest, ParentSetsTurnDictionaryEncodingOffForSetIdOnlyAndUseLargePages) {
    constexpr std::size_t n = 200'000;// 1.6 MB of set_id, more than one page at the writer's default 1 MiB
    Table parentSets{.columns = {{"iteration", I(n, 10)},
                                 {"phase", S(n, "sampling")},
                                 {"window", I(n, 1)},
                                 {"infection", S(n, "inf_a")},
                                 {"set_id", I(n)},
                                 {"prob", D(n, 0.5)}}};
    auto& setId = std::get<I>(parentSets.columns[4].second);
    for (std::size_t k = 0; k < n; ++k) {
        setId[k] = static_cast<std::int64_t>(k % 1000);
    }
    {
        ParquetStoreSink sink(root);
        sink.appendRows("parent_sets", parentSets);
        sink.appendRows("infections", infections(1, "sampling"));
        sink.commitWindow({.window = 1, .firstIteration = 1, .lastIteration = 10, .phase = "sampling"});
    }

    const auto reader = openParquet(root / "parent_sets" / "window-000001.parquet");
    const auto file = reader->parquet_reader();
    const auto meta = file->metadata();
    ASSERT_EQ(meta->num_row_groups(), 1);
    const auto rowGroup = meta->RowGroup(0);
    for (int c = 0; c < rowGroup->num_columns(); ++c) {
        const auto chunk = rowGroup->ColumnChunk(c);
        const auto name = meta->schema()->Column(c)->name();
        EXPECT_EQ(chunk->compression(), arrow::Compression::ZSTD) << name;
        EXPECT_EQ(chunk->has_dictionary_page(), name != "set_id") << name;
    }

    const int setIdColumn = meta->schema()->ColumnIndex("set_id");
    auto pages = file->RowGroup(0)->GetColumnPageReader(setIdColumn);
    int dataPages = 0;
    while (const auto page = pages->NextPage()) {
        dataPages += page->type() == parquet::PageType::DATA_PAGE or page->type() == parquet::PageType::DATA_PAGE_V2;
    }
    EXPECT_EQ(dataPages, 1);
    std::shared_ptr<arrow::Table> read;
    ASSERT_TRUE(reader->ReadTable(&read).ok());
    EXPECT_EQ(read->num_rows(), static_cast<std::int64_t>(n));

    const auto other = openParquet(root / "infections" / "window-000001.parquet")->parquet_reader()->metadata();
    const auto infection = other->RowGroup(0)->ColumnChunk(other->schema()->ColumnIndex("infection"));
    EXPECT_TRUE(infection->has_dictionary_page());
    EXPECT_EQ(infection->compression(), arrow::Compression::ZSTD);
}

// #55: `parent_sets.set_id` is delta-encoded. An Infection's ids come sorted within an iteration, so deltas are small.
TEST_F(ParquetStoreSinkTest, ParentSetsDeltaEncodeSetId) {
    constexpr std::size_t n = 1000;
    Table parentSets{.columns = {{"iteration", I(n, 10)},
                                 {"phase", S(n, "sampling")},
                                 {"window", I(n, 1)},
                                 {"infection", S(n, "inf_a")},
                                 {"set_id", I(n)},
                                 {"prob", D(n, 0.5)}}};
    auto& setId = std::get<I>(parentSets.columns[4].second);
    for (std::size_t k = 0; k < n; ++k) {
        setId[k] = static_cast<std::int64_t>(k);
    }
    {
        ParquetStoreSink sink(root);
        sink.appendRows("parent_sets", parentSets);
        sink.commitWindow({.window = 1, .firstIteration = 1, .lastIteration = 10, .phase = "sampling"});
    }

    const auto reader = openParquet(root / "parent_sets" / "window-000001.parquet");
    const auto meta   = reader->parquet_reader()->metadata();
    for (int g = 0; g < meta->num_row_groups(); ++g) {
        const auto chunk     = meta->RowGroup(g)->ColumnChunk(meta->schema()->ColumnIndex("set_id"));
        const auto encodings = chunk->encodings();
        EXPECT_NE(std::ranges::find(encodings, parquet::Encoding::DELTA_BINARY_PACKED), encodings.end());
        EXPECT_EQ(std::ranges::find(encodings, parquet::Encoding::PLAIN), encodings.end());
    }
    const auto read = readParquet(root / "parent_sets" / "window-000001.parquet");
    const auto ids  = std::static_pointer_cast<arrow::Int64Array>(read->GetColumnByName("set_id")->chunk(0));
    ASSERT_EQ(ids->length(), static_cast<std::int64_t>(n));
    for (std::size_t k = 0; k < n; ++k) {
        EXPECT_EQ(ids->Value(static_cast<std::int64_t>(k)), static_cast<std::int64_t>(k));
    }
}

TEST_F(ParquetStoreSinkDeathTest, UnwritableTableDirectoryExitsWithOutputFailureNamingThePath) {
    const auto file = root / "infections" / ".in-flight.parquet";
    EXPECT_EXIT(
            {
                ParquetStoreSink sink(root);
                sink.appendRows("infections", infections(0, "burn_in"));
                sink.commitWindow({.window = 0, .firstIteration = 0, .lastIteration = 10, .phase = "burn_in"});
                fs::permissions(root / "infections", fs::perms::owner_read | fs::perms::owner_exec);
                sink.appendRows("infections", infections(1, "sampling"));
            },
            ::testing::ExitedWithCode(transmission_nets::core::io::OUTPUT_FAILURE), "^ERROR: cannot open output file " + file.string() + ": ");
    fs::permissions(root / "infections", fs::perms::owner_all);
}

// The interrupt watcher may close a Window from its own thread while the main thread keeps running, so a failure
// there must end the process without running exit handlers the main thread could race.
TEST_F(ParquetStoreSinkDeathTest, UnwritableStoreFromAnotherThreadExitsWithoutExitHandlers) {
    fs::create_directories(root);
    fs::permissions(root, fs::perms::owner_read | fs::perms::owner_exec);
    EXPECT_EXIT(
            {
                std::atexit([] { std::fputs("exit handler ran\n", stderr); });
                ParquetStoreSink sink(root);
                std::thread([&] { sink.appendRows("infections", infections(0, "burn_in")); }).join();
                std::exit(0);
            },
            ::testing::ExitedWithCode(transmission_nets::core::io::OUTPUT_FAILURE),
            "^ERROR: cannot create the directory of output file " + (root / "infections").string() + ": [^\n]*\n$");
}

// A Run killed with a Window in flight reads back every Window it closed, and nothing of the one it didn't, though
// row groups of it already reached the disk.
TEST_F(ParquetStoreSinkDeathTest, KilledRunReadsBackEveryClosedWindowAndHidesTheInFlightOne) {
    EXPECT_EXIT(
            {
                ParquetStoreSink sink(root, 1);// every allowed end cuts a row group, so some of Window 2 is on disk
                sink.appendRows("infections", infections(0, "burn_in"));
                sink.commitWindow({.window = 0, .firstIteration = 0, .lastIteration = 10, .phase = "burn_in"});
                sink.appendRows("infections", infections(1, "sampling"));
                sink.commitWindow({.window = 1, .firstIteration = 1, .lastIteration = 10, .phase = "sampling"});
                sink.appendRows("infections", infections(2, "sampling"));
                sink.appendRows("chains", infections(2, "sampling"));
                sink.allowRowGroupEnd();
                sink.appendRows("infections", infections(2, "sampling"));
                std::raise(SIGKILL);
            },
            ::testing::KilledBySignal(SIGKILL), "");

    EXPECT_EQ(filesIn(root / "_closed"), (std::vector<std::string>{"window-000000", "window-000001"}));
    const auto closed = readClosed(root, "infections");
    ASSERT_EQ(closed.size(), 2);
    for (std::int64_t w = 0; w < 2; ++w) {
        ASSERT_EQ(closed[w]->num_rows(), 2);
        EXPECT_EQ(column<arrow::Int64Array>(*closed[w], "window")->Value(0), w);
    }
    for (const auto& name : filesIn(root / "infections")) {
        EXPECT_TRUE(name.starts_with(".") or name == "window-000000.parquet" or name == "window-000001.parquet") << name;
    }
    for (const auto& name : filesIn(root / "chains")) {
        EXPECT_TRUE(name.starts_with(".")) << name;
    }
    EXPECT_TRUE(readClosed(root, "chains").empty());
}
