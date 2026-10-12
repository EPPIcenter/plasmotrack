//
// Reads a Run's store back the way an analyst would: each table's rows over every closed Window, decoded into a Table.
// Shared by the CLI end-to-end tests.
//

#ifndef TRANSMISSION_NETWORKS_TEST_STOREREADING_H
#define TRANSMISSION_NETWORKS_TEST_STOREREADING_H

#include "core/io/Table.h"

#include <arrow/api.h>
#include <nlohmann/json.hpp>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace transmission_nets::test {
    namespace fs = std::filesystem;
    using core::io::Column;
    using core::io::Table;

    // One Parquet file, decoded column by column.
    inline Table readParquetTable(const fs::path& file) {
        auto in     = arrow::io::ReadableFile::Open(file.string()).ValueOrDie();
        auto reader = parquet::arrow::OpenFile(in, arrow::default_memory_pool()).ValueOrDie();
        std::shared_ptr<arrow::Table> table;
        if (!reader->ReadTable(&table).ok()) {
            throw std::runtime_error("cannot read " + file.string());
        }
        table = table->CombineChunks().ValueOrDie();

        Table out;
        for (int c = 0; c < table->num_columns(); ++c) {
            const auto name  = table->field(c)->name();
            const auto chunk = table->num_rows() == 0 ? nullptr : table->column(c)->chunk(0);
            const auto type  = table->field(c)->type()->id();
            const auto rows  = table->num_rows();
            Column column;
            if (type == arrow::Type::INT64) {
                std::vector<std::int64_t> v;
                for (std::int64_t r = 0; r < rows; ++r) v.push_back(std::static_pointer_cast<arrow::Int64Array>(chunk)->Value(r));
                column = std::move(v);
            } else if (type == arrow::Type::DOUBLE) {
                std::vector<double> v;
                for (std::int64_t r = 0; r < rows; ++r) v.push_back(std::static_pointer_cast<arrow::DoubleArray>(chunk)->Value(r));
                column = std::move(v);
            } else if (type == arrow::Type::STRING) {
                std::vector<std::string> v;
                for (std::int64_t r = 0; r < rows; ++r) v.push_back(std::static_pointer_cast<arrow::StringArray>(chunk)->GetString(r));
                column = std::move(v);
            } else if (type == arrow::Type::BOOL) {
                std::vector<bool> v;
                for (std::int64_t r = 0; r < rows; ++r) v.push_back(std::static_pointer_cast<arrow::BooleanArray>(chunk)->Value(r));
                column = std::move(v);
            } else {
                throw std::runtime_error("unexpected column type in " + file.string() + ": " + name);
            }
            out.columns.emplace_back(name, std::move(column));
        }
        return out;
    }

    // The Windows a reader may read: those with a commit marker, in order.
    inline std::vector<std::string> closedWindows(const fs::path& store) {
        std::vector<std::string> windows;
        if (fs::exists(store / "_closed")) {
            for (const auto& marker : fs::directory_iterator(store / "_closed")) {
                if (const auto name = marker.path().filename().string(); !name.starts_with(".")) {
                    windows.push_back(name);
                }
            }
        }
        std::ranges::sort(windows);
        return windows;
    }

    // A table's rows over every closed Window, in Window order. Empty if no closed Window holds the table.
    inline Table readStoreTable(const fs::path& store, const std::string& table) {
        Table out;
        for (const auto& window : closedWindows(store)) {
            if (const auto file = store / table / (window + ".parquet"); fs::exists(file)) {
                out.append(readParquetTable(file));
            }
        }
        return out;
    }

    // Every table a store holds: its Window tables over every closed Window, and its structure tables under
    // "structure/<name>".
    inline std::map<std::string, Table> readStore(const fs::path& store) {
        std::map<std::string, Table> out;
        for (const auto& entry : fs::directory_iterator(store)) {
            const auto name = entry.path().filename().string();
            if (!entry.is_directory() or name == "_closed" or name.starts_with(".")) {
                continue;
            }
            if (name == "structure") {
                for (const auto& file : fs::directory_iterator(entry.path())) {
                    if (file.path().extension() == ".parquet" and !file.path().filename().string().starts_with(".")) {
                        out["structure/" + file.path().stem().string()] = readParquetTable(file.path());
                    }
                }
            } else {
                out[name] = readStoreTable(store, name);
            }
        }
        return out;
    }

    // The events table without wall times or durations, which differ between Runs of one recipe.
    inline Table withoutWallTimes(const Table& events) {
        Table out;
        for (const auto& [name, column] : events.columns) {
            if (name == "wall_time") {
                continue;
            }
            if (name != "payload") {
                out.columns.emplace_back(name, column);
                continue;
            }
            std::vector<std::string> payloads;
            for (const auto& text : std::get<std::vector<std::string>>(column)) {
                auto payload = nlohmann::json::parse(text);
                for (const auto* key : {"started_at", "seconds", "iterations_per_second"}) {
                    payload.erase(key);
                }
                payloads.push_back(payload.dump());
            }
            out.columns.emplace_back(name, std::move(payloads));
        }
        return out;
    }
}// namespace transmission_nets::test

#endif//TRANSMISSION_NETWORKS_TEST_STOREREADING_H
