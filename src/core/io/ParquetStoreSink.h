//
// A store sink writing a Run's store as Parquet files under one directory:
//
//   <root>/manifest.json                    the Run manifest
//   <root>/structure/<table>.parquet        lookup tables
//   <root>/<table>/window-NNNNNN.parquet    one file per table per Window, none for a table with no rows
//   <root>/_closed/window-NNNNNN            the Window's commit marker, a JSON object
//
// Rows handed over for the in-flight Window wait in memory until their table's row group ends, then go to the table's
// file under a dot-prefixed name. A row group ends at the first point the caller allows once its rows take the
// row-group budget in memory, so a table that never reaches it writes one row group per Window. Committing the Window closes each file and moves it into place durably (sync, rename, sync the directory; see DurableFile.h), then writes the marker the same way. A reader counts a Window only once its marker exists, so a killed Run reads back up to its last closed
// Window.
//
// Files are zstd level 3. `parent_sets` writes `set_id` without a dictionary, as DELTA_BINARY_PACKED, and uses 8 MiB data
// pages (Arrow 22 has no rows-per-page setting); every other column keeps the writer's defaults.
//
// Any output failure ends the process with OUTPUT_FAILURE and a stderr message naming the path; see OutputFailure.h.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARQUETSTORESINK_H
#define TRANSMISSION_NETWORKS_APP_PARQUETSTORESINK_H

#include "OutputFailure.h"
#include "StoreSink.h"

#include <filesystem>
#include <map>
#include <memory>
#include <string>

namespace transmission_nets::core::io {
    class ParquetStoreSink final : public StoreSink {
    public:
        // The rows the sink holds per table, in bytes of Arrow memory, before it writes a row group. On the He01 Run of
        // #56 that is about 5 Sampling iterations of `parent_sets` and 17 of each genotype table, and the store is 2.5%
        // larger than with one row group per Window; 16 MiB makes it 11% larger, 8 MiB 25%.
        static constexpr std::int64_t ROW_GROUP_BYTES = std::int64_t{32} << 20;

        // Creates `root` if it does not exist.
        explicit ParquetStoreSink(std::filesystem::path root, std::int64_t rowGroupBytes = ROW_GROUP_BYTES);
        ~ParquetStoreSink() override;

        void writeManifest(const nlohmann::json& manifest) override;
        void writeStructure(const std::string& table, const Table& rows) override;
        void appendRows(const std::string& table, const Table& rows) override;
        void allowRowGroupEnd() override;
        void commitWindow(const WindowMarker& marker) override;

        // Removes markers first, so no reader sees a Window whose files are gone.
        void discardWindowsAfter(std::int64_t lastKept) override;

    private:
        [[nodiscard]] std::filesystem::path inFlightPath(const std::string& table) const;

        // A table's in-flight file, open until the Window commits.
        struct InFlightFile;

        std::filesystem::path root_;
        std::int64_t rowGroupBytes_;
        std::map<std::string, std::unique_ptr<InFlightFile>> inFlight_;
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_PARQUETSTORESINK_H
