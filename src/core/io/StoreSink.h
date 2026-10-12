//
// Where the Run recorder's store goes. The recorder decides which rows exist and which Window they belong to; the sink
// owns how tables are encoded and how a Window is committed.
//

#ifndef TRANSMISSION_NETWORKS_APP_STORESINK_H
#define TRANSMISSION_NETWORKS_APP_STORESINK_H

#include "Table.h"

#include <nlohmann/json.hpp>

#include <cstdint>
#include <string>

namespace transmission_nets::core::io {
    // What a Window's commit marker records.
    struct WindowMarker {
        std::int64_t window;
        std::int64_t firstIteration;
        std::int64_t lastIteration;
        std::string phase;

        bool operator==(const WindowMarker&) const = default;
    };

    class StoreSink {
    public:
        virtual ~StoreSink() = default;

        // The Run manifest, written once at Run start.
        virtual void writeManifest(const nlohmann::json& manifest) = 0;

        // A lookup table under `structure/`, written once at Run start.
        virtual void writeStructure(const std::string& table, const Table& rows) = 0;

        // Rows of a table for the in-flight Window, following those handed over before. Every call for a table has the
        // same columns. A call with no rows is ignored, and a table with none in the Window has nothing to commit.
        virtual void appendRows(const std::string& table, const Table& rows) = 0;

        // A row group may end here, after the rows handed over so far. A sink ends row groups only at these points and
        // at commit, and picks which of these it takes; row groups are how a table's file is cut, not what it holds.
        virtual void allowRowGroupEnd() = 0;

        // Makes the rows handed over since the last commit durable as one Window, then marks it closed.
        virtual void commitWindow(const WindowMarker& marker) = 0;

        // Before a Run resumes from a checkpoint taken once Window `lastKept` was committed: removes every later Window
        // and anything left in flight, so the store holds exactly what the checkpoint covers.
        virtual void discardWindowsAfter(std::int64_t lastKept) = 0;
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_STORESINK_H
