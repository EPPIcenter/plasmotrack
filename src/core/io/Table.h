//
// Rows of one store table, held column by column. A sink decides how a Table is encoded; the Run recorder only fills it.
//

#ifndef TRANSMISSION_NETWORKS_APP_TABLE_H
#define TRANSMISSION_NETWORKS_APP_TABLE_H

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace transmission_nets::core::io {
    using Column = std::variant<std::vector<std::int64_t>, std::vector<double>, std::vector<std::string>, std::vector<bool>>;

    struct Table {
        // Columns in order. Every column holds the same number of rows.
        std::vector<std::pair<std::string, Column>> columns{};

        [[nodiscard]] std::size_t rows() const {
            return columns.empty() ? 0 : std::visit([](const auto& c) { return c.size(); }, columns.front().second);
        }

        [[nodiscard]] std::vector<std::string> names() const {
            std::vector<std::string> out;
            for (const auto& [name, _] : columns) {
                out.push_back(name);
            }
            return out;
        }

        // The column called `name`, which must hold values of type T.
        template<typename T>
        [[nodiscard]] const std::vector<T>& column(const std::string& name) const {
            for (const auto& [n, c] : columns) {
                if (n == name) {
                    return std::get<std::vector<T>>(c);
                }
            }
            throw std::out_of_range("no column " + name);
        }

        // Appends the rows of `other`, which must have the same columns.
        void append(const Table& other) {
            if (columns.empty()) {
                columns = other.columns;
                return;
            }
            for (std::size_t k = 0; k < columns.size(); ++k) {
                std::visit([&](auto& to) {
                    const auto& from = std::get<std::decay_t<decltype(to)>>(other.columns.at(k).second);
                    to.insert(to.end(), from.begin(), from.end());
                },
                           columns[k].second);
            }
        }
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_TABLE_H
