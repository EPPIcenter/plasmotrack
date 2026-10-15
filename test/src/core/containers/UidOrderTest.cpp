//
// Parent sets and per-locus maps must iterate in construction (uid) order, not address order, so that a run is
// reproducible for a given seed.
//

#include "gtest/gtest.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/ParentSet.h"
#include "core/datatypes/Genotype.h"

#include <array>
#include <memory>
#include <new>
#include <string>
#include <vector>

using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;

namespace {
    // Constructs objects into slots from last to first, so later-constructed objects have lower addresses.
    template<typename T, size_t N>
    class DescendingAddresses {
    public:
        template<typename Make>
        explicit DescendingAddresses(Make&& make) {
            for (size_t i = 0; i < N; ++i) {
                void* slot = &storage_[N - 1 - i];
                objects_.push_back(std::shared_ptr<T>(make(slot, i), [](T* object) { object->~T(); }));
            }
        }
        const std::vector<std::shared_ptr<T>>& inConstructionOrder() const { return objects_; }

    private:
        struct alignas(T) Slot {
            std::byte bytes[sizeof(T)];
        };
        std::array<Slot, N> storage_{};
        std::vector<std::shared_ptr<T>> objects_{};
    };
}// namespace

TEST(UidOrderTest, ParentSetsIterateInConstructionOrder) {
    using InfectionEvent = Infection<>;
    constexpr size_t n   = 6;
    DescendingAddresses<InfectionEvent, n> infections([](void* slot, const size_t i) {
        return new (slot) InfectionEvent("inf" + std::to_string(i), 10.0 + i);
    });
    const auto& expected = infections.inConstructionOrder();
    ASSERT_GT(expected.front().get(), expected.back().get());

    ParentSet<InfectionEvent> parents{};
    for (const size_t i : {3, 0, 5, 1, 4, 2}) {
        parents.insert(expected[i]);
    }
    ASSERT_EQ(parents.size(), n);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(parents.begin()[i], expected[i]) << "position " << i;
    }
}

TEST(UidOrderTest, GenotypeMapsIterateInConstructionOrder) {
    using InfectionEvent = Infection<>;
    constexpr size_t n   = 5;
    DescendingAddresses<Locus, n> loci([](void* slot, const size_t i) {
        return new (slot) Locus("L" + std::to_string(i), 4);
    });
    const auto& expected = loci.inConstructionOrder();
    ASSERT_GT(expected.front().get(), expected.back().get());

    InfectionEvent infection("inf", 10.0);
    for (const size_t i : {2, 4, 0, 3, 1}) {
        infection.addLatentGenetics(expected[i], Genotype("0100"));
    }
    // loci() keeps insertion order; the per-locus genotype maps are ordered by key.
    std::vector<std::shared_ptr<Locus>> order;
    for (const auto& [locus, genotype] : infection.latentGenotype()) {
        order.push_back(locus);
    }
    ASSERT_EQ(order.size(), n);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(order[i], expected[i]) << "position " << i;
    }
}
