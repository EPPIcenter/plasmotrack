//
// Orders objects by their uid instead of their address.
//

#ifndef TRANSMISSION_NETWORKS_APP_UIDLESS_H
#define TRANSMISSION_NETWORKS_APP_UIDLESS_H

#include <functional>
#include <memory>

namespace transmission_nets::core::containers {

    /*
     * Containers that iterate in key order (parent sets, per-locus maps) decide the order in which parent sets are
     * enumerated and sampled and in which per-locus terms are summed. Ordering by address makes that order, and with
     * it the whole chain, differ from run to run for the same seed; uids are assigned in construction order, which
     * follows the input.
     */
    struct UidLess {
        // Types without a uid (e.g. in tests) fall back to address order.
        template<typename T>
        bool operator()(const T* lhs, const T* rhs) const noexcept {
            if constexpr (requires { lhs->uid(); }) {
                return lhs->uid() < rhs->uid();
            } else if constexpr (requires { lhs->uid; }) {
                return lhs->uid < rhs->uid;
            } else {
                return std::less<const T*>{}(lhs, rhs);
            }
        }

        template<typename T>
        bool operator()(const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const noexcept {
            return (*this)(lhs.get(), rhs.get());
        }
    };

}// namespace transmission_nets::core::containers

#endif//TRANSMISSION_NETWORKS_APP_UIDLESS_H
