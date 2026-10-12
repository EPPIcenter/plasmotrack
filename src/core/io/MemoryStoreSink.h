//
// A store sink that keeps everything in memory, for tests. Rows handed over for a Window are visible only once the
// Window is committed, as a reader of a store on disk would see them.
//

#ifndef TRANSMISSION_NETWORKS_APP_MEMORYSTORESINK_H
#define TRANSMISSION_NETWORKS_APP_MEMORYSTORESINK_H

#include "StoreSink.h"

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace transmission_nets::core::io {
    class MemoryStoreSink final : public StoreSink {
    public:
        void writeManifest(const nlohmann::json& manifest) override { manifest_ = manifest; }

        void writeStructure(const std::string& table, const Table& rows) override { structure_[table] = rows; }

        void appendRows(const std::string& table, const Table& rows) override {
            if (rows.rows() == 0) {
                return;
            }
            auto& inFlight = inFlight_[table];
            if (!inFlight.columns.empty() and inFlight.names() != rows.names()) {
                throw std::logic_error(table + " handed over with different columns in one Window");
            }
            inFlight.append(rows);
            pending_[table] += rows.rows();
        }

        // Ends a row group at every point allowed.
        void allowRowGroupEnd() override {
            for (const auto& [table, rows] : pending_) {
                inFlightGroups_[table].push_back(rows);
            }
            pending_.clear();
        }

        void commitWindow(const WindowMarker& marker) override {
            allowRowGroupEnd();
            committed_.push_back({.marker = marker, .tables = std::move(inFlight_), .rowGroups = std::move(inFlightGroups_)});
            inFlight_.clear();
            inFlightGroups_.clear();
        }

        void discardWindowsAfter(const std::int64_t lastKept) override {
            std::erase_if(committed_, [&](const Committed& w) { return w.marker.window > lastKept; });
            inFlight_.clear();
            inFlightGroups_.clear();
            pending_.clear();
        }

        [[nodiscard]] const nlohmann::json& manifest() const { return manifest_; }

        [[nodiscard]] const Table& structure(const std::string& table) const { return structure_.at(table); }

        // The markers of every committed Window, in commit order.
        [[nodiscard]] std::vector<WindowMarker> windows() const {
            std::vector<WindowMarker> out;
            for (const auto& w : committed_) {
                out.push_back(w.marker);
            }
            return out;
        }

        // The tables handed over for one committed Window, by its position in commit order.
        [[nodiscard]] const std::map<std::string, Table>& windowTables(const std::size_t k) const { return committed_.at(k).tables; }

        // The rows in each row group of a table in one committed Window, in the order handed over.
        [[nodiscard]] const std::vector<std::size_t>& rowGroups(const std::size_t k, const std::string& table) const { return committed_.at(k).rowGroups.at(table); }

        // A table's rows across every committed Window, in commit order.
        [[nodiscard]] Table table(const std::string& name) const {
            Table out;
            for (const auto& w : committed_) {
                if (const auto found = w.tables.find(name); found != w.tables.end()) {
                    out.append(found->second);
                }
            }
            return out;
        }

        // Rows handed over for a Window not yet committed.
        [[nodiscard]] bool hasUncommittedTables() const { return !inFlight_.empty(); }

        // A table's rows handed over for a Window not yet committed.
        [[nodiscard]] std::size_t uncommittedRows(const std::string& table) const {
            const auto found = inFlight_.find(table);
            return found == inFlight_.end() ? 0 : found->second.rows();
        }

    private:
        struct Committed {
            WindowMarker marker;
            std::map<std::string, Table> tables;
            std::map<std::string, std::vector<std::size_t>> rowGroups;
        };

        nlohmann::json manifest_{};
        std::map<std::string, Table> structure_{};
        std::map<std::string, Table> inFlight_{};
        std::map<std::string, std::vector<std::size_t>> inFlightGroups_{};
        std::map<std::string, std::size_t> pending_{};// rows per table since the last row group ended
        std::vector<Committed> committed_{};
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_MEMORYSTORESINK_H
