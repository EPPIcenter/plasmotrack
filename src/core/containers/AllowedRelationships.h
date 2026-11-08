//
// Created by Maxwell Murphy on 1/10/24.
//

#ifndef ALLOWEDRELATIONSHIPS_H
#define ALLOWEDRELATIONSHIPS_H

#include <memory>
#include <vector>

#include <boost/container/flat_map.hpp>

#include "core/containers/UidLess.h"

namespace transmission_nets::core::containers {
    template<typename InfectionEvent>
    struct AllowedRelationships {
        const std::vector<std::shared_ptr<InfectionEvent>>& permittedParents(const std::shared_ptr<InfectionEvent>& infectionEvent) const {
            const auto it = permittedParents_.find(infectionEvent);
            if (it == permittedParents_.end()) {
                return empty_;
            }
            return it->second;
        }

        const std::vector<std::shared_ptr<InfectionEvent>>& allowedChildren(const std::shared_ptr<InfectionEvent>& infectionEvent) const {
            const auto it = allowedChildren_.find(infectionEvent);
            if (it == allowedChildren_.end()) {
                return empty_;
            }
            return it->second;
        }

        void addParent(const std::shared_ptr<InfectionEvent>& infectionEvent, const std::shared_ptr<InfectionEvent>& parent) {
            permittedParents_[infectionEvent].push_back(parent);
        }

        void addChild(const std::shared_ptr<InfectionEvent>& infectionEvent, const std::shared_ptr<InfectionEvent>& child) {
            allowedChildren_[infectionEvent].push_back(child);
        }

    private:
        inline static const std::vector<std::shared_ptr<InfectionEvent>> empty_{};
        boost::container::flat_map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>, UidLess> permittedParents_{};
        boost::container::flat_map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>, UidLess> allowedChildren_{};
    };
}// namespace transmission_nets::core::containers

#endif//ALLOWEDRELATIONSHIPS_H
