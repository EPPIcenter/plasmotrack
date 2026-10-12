#ifndef TRANSMISSION_NETWORKS_APP_SAMPLERINDEX_H
#define TRANSMISSION_NETWORKS_APP_SAMPLERINDEX_H

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/samplers/AbstractSampler.h"

namespace transmission_nets::core::samplers {

    // The samplers registered with one scheduler, by id. Ids are `<type>/<target>` and unique within the scheduler.
    class SamplerIndex {
    public:
        // Throws std::invalid_argument when `id` is already registered.
        void add(const std::string& id, AbstractSampler& sampler) {
            if (!by_id_.emplace(id, &sampler).second) {
                throw std::invalid_argument("sampler id registered twice: " + id);
            }
            ids_.push_back(id);
        }

        // In registration order.
        [[nodiscard]] const std::vector<std::string>& ids() const noexcept { return ids_; }

        // Throws std::out_of_range for an id that isn't registered.
        [[nodiscard]] AbstractSampler& at(const std::string& id) const {
            const auto it = by_id_.find(id);
            if (it == by_id_.end()) {
                throw std::out_of_range("no sampler with id " + id);
            }
            return *it->second;
        }

    private:
        std::vector<std::string> ids_{};
        std::unordered_map<std::string, AbstractSampler*> by_id_{};
    };

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERINDEX_H
