//
// Created by Maxwell Murphy on 1/10/24.
//

#ifndef ALLOWEDRELATIONSHIPS_H
#define ALLOWEDRELATIONSHIPS_H

// #include <map>
#include <memory>
#include <ranges>
#include <vector>

#include <fmt/format.h>

#include <boost/container/flat_map.hpp>


namespace transmission_nets::core::containers {
    template<typename InfectionEvent>
    struct AllowedRelationships {
        std::vector<std::shared_ptr<InfectionEvent>> allowedParents(const std::shared_ptr<InfectionEvent> infectionEvent) {
            if (allowedParents_.find(infectionEvent) == allowedParents_.end()) {
                return {};
            }
            std::vector<std::shared_ptr<InfectionEvent>> out;
            for (const auto& parent : allowedParents_.at(infectionEvent)) {
                out.push_back(infectionEvents_[parent]);
            }
            return out;
        }

        std::vector<std::shared_ptr<InfectionEvent>> allowedChildren(const std::shared_ptr<InfectionEvent> infectionEvent) {
            if (allowedChildren_.find(infectionEvent) == allowedChildren_.end()) {
                return {};
            }
            std::vector<std::shared_ptr<InfectionEvent>> out;
            for (const auto& child : allowedChildren_.at(infectionEvent)) {
                out.push_back(infectionEvents_[child]);
            }
            return out;
        }

        void addParent(const std::shared_ptr<InfectionEvent> infectionEvent, const std::shared_ptr<InfectionEvent> parent) {
            if (allowedParents_.find(infectionEvent) == allowedParents_.end()) {
                allowedParents_.emplace(infectionEvent, std::vector<unsigned short>{});
            }
            allowedParents_.at(infectionEvent).push_back(parent->uid());
        }

        void addChild(const std::shared_ptr<InfectionEvent> infectionEvent, const std::shared_ptr<InfectionEvent> child) {
            if (allowedChildren_.find(infectionEvent) == allowedChildren_.end()) {
                allowedChildren_.emplace(infectionEvent, std::vector<unsigned short>{});
            }
            allowedChildren_.at(infectionEvent).push_back(child->uid());
        }

        void addInfectionEvent(const std::shared_ptr<InfectionEvent> infectionEvent) {
            if (static_cast<unsigned short>(infectionEvents_.size()) < infectionEvent->uid() + 1) {
                infectionEvents_.resize(infectionEvent->uid() + 1);
            }
            infectionEvents_[infectionEvent->uid()] = infectionEvent;
        }

    private:
        std::vector<std::shared_ptr<InfectionEvent>> infectionEvents_{};
        boost::container::flat_map<std::shared_ptr<InfectionEvent>, std::vector<unsigned short>> allowedParents_{};
        boost::container::flat_map<std::shared_ptr<InfectionEvent>, std::vector<unsigned short>> allowedChildren_{};
    };
}// namespace transmission_nets::core::containers


#endif//ALLOWEDRELATIONSHIPS_H
