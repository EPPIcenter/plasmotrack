//
// The Parquet store sink. See ParquetStoreSink.h.
//

#include "ParquetStoreSink.h"
#include "DurableFile.h"
#include "OutputFailure.h"

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <fmt/core.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>


#include <algorithm>
#include <charconv>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <system_error>

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;

    namespace {
        // Pages of 2^20 rows for the 8-byte columns of `parent_sets`. Arrow C++ cuts pages by bytes, not rows.
        constexpr std::int64_t PARENT_SETS_PAGE_BYTES = (std::int64_t{1} << 20) * 8;

        // At least six digits; a larger number takes as many as it needs.
        std::string windowName(const std::int64_t window) {
            return fmt::format("window-{:06d}", window);
        }

        // The Window a file or marker named by windowName belongs to: the digits after `window-`, of any width, up to
        // the extension if there is one. None for any other name.
        std::optional<std::int64_t> windowNumber(const std::string_view name) {
            constexpr std::string_view prefix = "window-";
            if (!name.starts_with(prefix)) {
                return std::nullopt;
            }
            auto digits = name.substr(prefix.size());
            digits      = digits.substr(0, digits.find('.'));
            if (digits.empty() or !std::ranges::all_of(digits, [](const char c) { return c >= '0' and c <= '9'; })) {
                return std::nullopt;
            }
            std::int64_t window  = 0;
            const auto [end, ec] = std::from_chars(digits.data(), digits.data() + digits.size(), window);
            if (ec != std::errc{} or end != digits.data() + digits.size()) {
                return std::nullopt;
            }
            return window;
        }

        void check(const arrow::Status& status, const fs::path& path) {
            if (!status.ok()) {
                failOutput("write", path, status.ToString());
            }
        }

        // An empty Arrow builder for a Table column's type.
        struct NewBuilder {
            std::unique_ptr<arrow::ArrayBuilder> operator()(const std::vector<std::int64_t>&) const { return std::make_unique<arrow::Int64Builder>(); }
            std::unique_ptr<arrow::ArrayBuilder> operator()(const std::vector<double>&) const { return std::make_unique<arrow::DoubleBuilder>(); }
            std::unique_ptr<arrow::ArrayBuilder> operator()(const std::vector<std::string>&) const { return std::make_unique<arrow::StringBuilder>(); }
            std::unique_ptr<arrow::ArrayBuilder> operator()(const std::vector<bool>&) const { return std::make_unique<arrow::BooleanBuilder>(); }
        };

        // Appends a Table column's values to the builder NewBuilder made for it.
        struct AppendTo {
            arrow::ArrayBuilder& builder;

            arrow::Status operator()(const std::vector<std::int64_t>& v) const { return static_cast<arrow::Int64Builder&>(builder).AppendValues(v); }
            arrow::Status operator()(const std::vector<double>& v) const { return static_cast<arrow::DoubleBuilder&>(builder).AppendValues(v); }
            arrow::Status operator()(const std::vector<std::string>& v) const { return static_cast<arrow::StringBuilder&>(builder).AppendValues(v); }
            arrow::Status operator()(const std::vector<bool>& v) const { return static_cast<arrow::BooleanBuilder&>(builder).AppendValues(v); }
        };

        // What a Table column's values take in an Arrow array.
        struct ArrowBytes {
            template<typename T>
            std::int64_t operator()(const std::vector<T>& v) const { return static_cast<std::int64_t>(v.size() * sizeof(T)); }
            std::int64_t operator()(const std::vector<bool>& v) const { return static_cast<std::int64_t>((v.size() + 7) / 8); }
            std::int64_t operator()(const std::vector<std::string>& v) const {
                std::int64_t out = 0;
                for (const auto& s : v) {
                    out += static_cast<std::int64_t>(sizeof(std::int32_t) + s.size());
                }
                return out;
            }
        };

        std::shared_ptr<parquet::WriterProperties> writerProperties(const std::string& table) {
            parquet::WriterProperties::Builder builder;
            builder.compression(arrow::Compression::ZSTD)->compression_level(3);
            if (table == "parent_sets") {
                // Sorted within an Infection (see RunRecorder), set ids delta-encode to a few bits each (#55).
                builder.disable_dictionary("set_id")->encoding("set_id", parquet::Encoding::DELTA_BINARY_PACKED)->data_pagesize(PARENT_SETS_PAGE_BYTES);
            }
            return builder.build();
        }
    }// namespace

    // A table's Parquet file open for writing. Rows wait in Arrow builders, more compact than a Table, until their row
    // group ends.
    struct ParquetStoreSink::InFlightFile {
        fs::path path;
        std::vector<std::unique_ptr<arrow::ArrayBuilder>> builders;
        std::shared_ptr<arrow::Schema> schema;
        std::shared_ptr<arrow::io::FileOutputStream> out;
        std::unique_ptr<parquet::arrow::FileWriter> writer;
        std::int64_t pending = 0;// rows in the builders
        std::int64_t bytes   = 0;// what they take in Arrow memory: values, string offsets and bytes

        // Opens `p` for rows with the columns of `like`.
        InFlightFile(const std::string& table, const Table& like, fs::path p) : path(std::move(p)) {
            arrow::FieldVector fields;
            for (const auto& [name, values] : like.columns) {
                builders.push_back(std::visit(NewBuilder{}, values));
                fields.push_back(arrow::field(name, builders.back()->type()));
            }
            schema = arrow::schema(fields);

            auto stream = arrow::io::FileOutputStream::Open(path.string());
            if (!stream.ok()) {
                failOutput("open", path, stream.status().ToString());
            }
            out           = *stream;
            auto parquet_ = parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), out, writerProperties(table));
            if (!parquet_.ok()) {
                failOutput("write", path, parquet_.status().ToString());
            }
            writer = std::move(*parquet_);
        }

        void append(const Table& rows) {
            if (rows.names() != schema->field_names()) {
                throw std::logic_error(path.string() + ": rows handed over with different columns in one Window");
            }
            for (std::size_t c = 0; c < builders.size(); ++c) {
                check(std::visit(AppendTo{*builders[c]}, rows.columns[c].second), path);
                bytes += std::visit(ArrowBytes{}, rows.columns[c].second);
            }
            pending += static_cast<std::int64_t>(rows.rows());
        }

        // Writes the rows in the builders as one row group. If `more` rows will follow, the builders then reserve room
        // for a group a little larger, so the next one doesn't grow them by copying.
        void endRowGroup(const bool more = true) {
            if (pending == 0) {
                return;
            }
            arrow::ArrayVector arrays;
            std::vector<std::int64_t> stringBytes;
            for (const auto& builder : builders) {
                const auto* strings = dynamic_cast<arrow::StringBuilder*>(builder.get());
                stringBytes.push_back(strings ? strings->value_data_length() : 0);
                std::shared_ptr<arrow::Array> array;
                check(builder->Finish(&array), path);
                arrays.push_back(std::move(array));
            }
            const auto rows = arrow::Table::Make(schema, std::move(arrays), pending);
            check(writer->WriteTable(*rows, pending), path);

            const auto room = [](const std::int64_t n) { return n + n / 8; };
            for (std::size_t c = 0; more and c < builders.size(); ++c) {
                check(builders[c]->Reserve(room(pending)), path);
                if (auto* strings = dynamic_cast<arrow::StringBuilder*>(builders[c].get())) {
                    check(strings->ReserveData(room(stringBytes[c])), path);
                }
            }
            pending = 0;
            bytes   = 0;
        }

        void close() {
            endRowGroup(false);
            check(writer->Close(), path);
            check(out->Close(), path);
        }
    };

    ParquetStoreSink::ParquetStoreSink(fs::path root, const std::int64_t rowGroupBytes) : root_(std::move(root)), rowGroupBytes_(rowGroupBytes) {
        createDirectories(root_);
    }

    // Files still in flight keep their dot-prefixed names; a resumed Run removes them.
    ParquetStoreSink::~ParquetStoreSink() = default;

    void ParquetStoreSink::writeManifest(const nlohmann::json& manifest) {
        writeDurably(root_ / "manifest.json", [&](const fs::path& p) { writeFile(p, manifest.dump(2) + "\n"); });
    }

    void ParquetStoreSink::writeStructure(const std::string& table, const Table& rows) {
        writeDurably(root_ / "structure" / (table + ".parquet"), [&](const fs::path& p) {
            InFlightFile file(table, rows, p);
            file.append(rows);
            file.close();
        });
    }

    void ParquetStoreSink::appendRows(const std::string& table, const Table& rows) {
        if (rows.rows() == 0) {
            return;
        }
        auto& file = inFlight_[table];
        if (!file) {
            createDirectories(root_ / table);
            file = std::make_unique<InFlightFile>(table, rows, inFlightPath(table));
        }
        file->append(rows);
    }

    void ParquetStoreSink::allowRowGroupEnd() {
        for (const auto& [_, file] : inFlight_) {
            if (file->bytes >= rowGroupBytes_) {
                file->endRowGroup();
            }
        }
    }

    void ParquetStoreSink::commitWindow(const WindowMarker& marker) {
        const auto name = windowName(marker.window);
        for (const auto& [table, file] : inFlight_) {
            file->close();
            moveDurably(file->path, root_ / table / (name + ".parquet"));
        }
        // A table's directory may be new in this Window.
        sync(root_);
        inFlight_.clear();

        const nlohmann::json closed{
                {"window", marker.window},
                {"first_iteration", marker.firstIteration},
                {"last_iteration", marker.lastIteration},
                {"phase", marker.phase},
        };
        writeDurably(root_ / "_closed" / name, [&](const fs::path& p) { writeFile(p, closed.dump() + "\n"); });
    }

    void ParquetStoreSink::discardWindowsAfter(const std::int64_t lastKept) {
        // Hidden files are in flight: a table not yet committed, or a marker not yet in place.
        const auto discarded = [&](const fs::path& file) {
            const auto name = file.filename().string();
            if (name.starts_with(".")) {
                return true;
            }
            const auto window = windowNumber(name);
            return window and *window > lastKept;
        };
        const auto discardIn = [&](const fs::path& dir) {
            if (!fs::is_directory(dir)) {
                return;
            }
            for (const auto& entry : fs::directory_iterator(dir)) {
                if (discarded(entry.path())) {
                    std::error_code ec;
                    fs::remove(entry.path(), ec);
                    if (ec) {
                        failOutput("remove", entry.path(), ec.message());
                    }
                }
            }
            sync(dir);
        };
        // Markers first, so no reader sees a Window whose files are gone.
        discardIn(root_ / "_closed");
        for (const auto& entry : fs::directory_iterator(root_)) {
            const auto name = entry.path().filename().string();
            if (entry.is_directory() and name != "_closed" and name != "structure") {
                discardIn(entry.path());
            }
        }
    }

    fs::path ParquetStoreSink::inFlightPath(const std::string& table) const {
        return root_ / table / ".in-flight.parquet";
    }
}// namespace transmission_nets::core::io
