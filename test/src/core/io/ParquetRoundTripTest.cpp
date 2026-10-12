//
// The build links Arrow C++ with Parquet and zstd: a small table written as a zstd-compressed Parquet file reads
// back with the same schema and values.
//

#include "gtest/gtest.h"

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <unistd.h>

#include <filesystem>
#include <memory>
#include <string>

namespace fs = std::filesystem;

namespace {

std::shared_ptr<arrow::Table> makeTable() {
    arrow::Int64Builder iteration;
    arrow::StringBuilder infection;
    arrow::DoubleBuilder prob;
    EXPECT_TRUE(iteration.AppendValues({0, 5, 10}).ok());
    EXPECT_TRUE(infection.AppendValues({"inf_a", "inf_b", "inf_a"}).ok());
    EXPECT_TRUE(prob.AppendValues({1.0, 0.25, 1e-9}).ok());

    auto schema = arrow::schema({arrow::field("iteration", arrow::int64()),
                                 arrow::field("infection", arrow::utf8()),
                                 arrow::field("prob", arrow::float64())});
    return arrow::Table::Make(schema, {iteration.Finish().ValueOrDie(), infection.Finish().ValueOrDie(),
                                       prob.Finish().ValueOrDie()});
}

}// namespace

TEST(ParquetRoundTripTest, ZstdParquetFileReadsBackTheTableItWasWrittenFrom) {
    const fs::path dir = fs::temp_directory_path() / ("plasmotrack-parquet-" + std::to_string(::getpid()));
    fs::remove_all(dir);
    fs::create_directories(dir);
    const auto file = (dir / "window-000001.parquet").string();

    const auto written = makeTable();
    {
        auto out = arrow::io::FileOutputStream::Open(file).ValueOrDie();
        auto props = parquet::WriterProperties::Builder()
                             .compression(arrow::Compression::ZSTD)
                             ->compression_level(3)
                             ->build();
        ASSERT_TRUE(parquet::arrow::WriteTable(*written, arrow::default_memory_pool(), out, 1024, props).ok());
        ASSERT_TRUE(out->Close().ok());
    }

    auto in = arrow::io::ReadableFile::Open(file).ValueOrDie();
    auto reader = parquet::arrow::OpenFile(in, arrow::default_memory_pool()).ValueOrDie();
    EXPECT_EQ(reader->parquet_reader()->metadata()->RowGroup(0)->ColumnChunk(0)->compression(),
              arrow::Compression::ZSTD);
    std::shared_ptr<arrow::Table> read;
    ASSERT_TRUE(reader->ReadTable(&read).ok());

    EXPECT_TRUE(read->schema()->Equals(*written->schema()));
    EXPECT_EQ(read->num_rows(), 3);
    auto prob = std::static_pointer_cast<arrow::DoubleArray>(read->GetColumnByName("prob")->chunk(0));
    EXPECT_DOUBLE_EQ(prob->Value(1), 0.25);
    EXPECT_DOUBLE_EQ(prob->Value(2), 1e-9);
    auto infection = std::static_pointer_cast<arrow::StringArray>(read->GetColumnByName("infection")->chunk(0));
    EXPECT_EQ(infection->GetString(1), "inf_b");

    fs::remove_all(dir);
}
