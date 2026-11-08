//
// A value shared between copies until one of them writes it.
//

#ifndef TRANSMISSION_NETWORKS_APP_COPYONWRITE_H
#define TRANSMISSION_NETWORKS_APP_COPYONWRITE_H

#include <cstdint>

namespace transmission_nets::core::containers {

    // Copying copies a pointer; mutate() clones the value first if another copy still holds it. Lets a saved state's
    // copy of a cache share every value it does not change. Not thread-safe (a chain's terms are used from one thread).
    template<typename V>
    class CopyOnWrite {
        struct Node {
            V value;
            std::uint32_t refs;
        };
        Node* node_ = nullptr;

    public:
        CopyOnWrite() = default;
        CopyOnWrite(const CopyOnWrite& other) noexcept : node_(other.node_) {
            if (node_ != nullptr) {
                ++node_->refs;
            }
        }
        CopyOnWrite(CopyOnWrite&& other) noexcept : node_(other.node_) { other.node_ = nullptr; }
        CopyOnWrite& operator=(const CopyOnWrite& other) noexcept {
            if (other.node_ != nullptr) {
                ++other.node_->refs;
            }
            reset();
            node_ = other.node_;
            return *this;
        }
        CopyOnWrite& operator=(CopyOnWrite&& other) noexcept {
            if (this != &other) {
                reset();
                node_       = other.node_;
                other.node_ = nullptr;
            }
            return *this;
        }
        ~CopyOnWrite() { reset(); }

        // The value; a default one if never set.
        [[nodiscard]] const V& get() const noexcept {
            static const V empty{};
            return node_ != nullptr ? node_->value : empty;
        }
        // The value, this copy's own.
        V& mutate() {
            if (node_ == nullptr) {
                node_ = new Node{V{}, 1};
            } else if (node_->refs > 1) {
                --node_->refs;
                node_ = new Node{node_->value, 1};
            }
            return node_->value;
        }
        void set(const V& value) {
            if (node_ != nullptr and node_->refs == 1) {
                node_->value = value;
            } else {
                reset();
                node_ = new Node{value, 1};
            }
        }
        // Lets go of the value, so that other copies holding it no longer share it with this one.
        void reset() noexcept {
            if (node_ != nullptr and --node_->refs == 0) {
                delete node_;
            }
            node_ = nullptr;
        }
    };

}// namespace transmission_nets::core::containers

#endif//TRANSMISSION_NETWORKS_APP_COPYONWRITE_H
