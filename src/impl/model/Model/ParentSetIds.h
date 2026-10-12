//
// The integer ids a Run gives its Parent sets: per Infection, numbered in order of first appearance. The Run recorder
// owns the dictionary while a Run goes on; a checkpoint carries it to a resumed Run.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETIDS_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETIDS_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {
    // One Infection's Parent-set ids: each Parent set, named by its members' ids in ascending order, and the id it was
    // given, plus the id the next new Parent set gets.
    struct ParentSetIds {
        std::map<std::vector<std::string>, std::uint64_t> ids;
        std::uint64_t nextId = 0;

        bool operator==(const ParentSetIds&) const = default;
    };

    // By Infection id.
    using ParentSetIdDictionary = std::map<std::string, ParentSetIds>;
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETIDS_H
