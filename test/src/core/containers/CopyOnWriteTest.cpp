//
// A CopyOnWrite value is shared between copies until one of them writes it, and a write never shows in another copy.
//

#include "gtest/gtest.h"

#include "core/containers/CopyOnWrite.h"

#include <array>
#include <utility>
#include <vector>

using transmission_nets::core::containers::CopyOnWrite;

namespace {
    using Value = std::array<double, 4>;

    CopyOnWrite<Value> holding(const double x) {
        CopyOnWrite<Value> box;
        box.set(Value{x, x, x, x});
        return box;
    }
}// namespace

TEST(CopyOnWriteTest, NeverSetReadsTheDefault) {
    const CopyOnWrite<Value> box;
    EXPECT_EQ(box.get(), Value{});
}

TEST(CopyOnWriteTest, CopiesShare) {
    const auto a = holding(1.0);
    const auto b = a;
    EXPECT_EQ(&a.get(), &b.get());
}

TEST(CopyOnWriteTest, MutatingASharedCopyClonesAndLeavesTheOtherIntact) {
    auto a           = holding(1.0);
    auto b           = a;
    const auto* held = &a.get();
    b.mutate()[0]    = 2.0;
    EXPECT_EQ(a.get(), (Value{1.0, 1.0, 1.0, 1.0}));
    EXPECT_EQ(b.get(), (Value{2.0, 1.0, 1.0, 1.0}));
    EXPECT_EQ(&a.get(), held);
    EXPECT_NE(&b.get(), held);
    // a now holds its value alone, so writing it does not clone.
    a.mutate()[1] = 3.0;
    EXPECT_EQ(&a.get(), held);
}

TEST(CopyOnWriteTest, MutatingAnOwnValueWritesInPlace) {
    auto a           = holding(1.0);
    const auto* held = &a.get();
    a.mutate()[0]    = 2.0;
    EXPECT_EQ(&a.get(), held);
    EXPECT_EQ(a.get()[0], 2.0);
}

TEST(CopyOnWriteTest, SetReleasesASharedValue) {
    auto a           = holding(1.0);
    auto b           = a;
    const auto* held = &a.get();
    b.set(Value{5.0, 5.0, 5.0, 5.0});
    EXPECT_EQ(a.get(), (Value{1.0, 1.0, 1.0, 1.0}));
    EXPECT_EQ(b.get(), (Value{5.0, 5.0, 5.0, 5.0}));
    // a is the value's only holder again.
    a.mutate()[0] = 2.0;
    EXPECT_EQ(&a.get(), held);
}

TEST(CopyOnWriteTest, ResetReleasesASharedValue) {
    auto a           = holding(1.0);
    auto b           = a;
    const auto* held = &a.get();
    b.reset();
    EXPECT_EQ(b.get(), Value{});
    a.mutate()[0] = 2.0;
    EXPECT_EQ(&a.get(), held);
}

// Saved states as the Transmission term keeps them: a save copies the state below, a change writes the copy, a
// restore drops back to the one below and an accept keeps the copy.
TEST(CopyOnWriteTest, SaveRestoreAndAcceptDoNotAlias) {
    std::vector<CopyOnWrite<Value>> states(2);
    states[0].set(Value{1.0, 1.0, 1.0, 1.0});

    // Save, change, restore: the saved state is untouched.
    states[1]          = states[0];
    states[1].mutate() = Value{2.0, 2.0, 2.0, 2.0};
    EXPECT_EQ(states[0].get(), (Value{1.0, 1.0, 1.0, 1.0}));

    // Save, change, accept: the change is kept, and the state given up shares nothing with it.
    states[1]          = states[0];
    states[1].mutate() = Value{3.0, 3.0, 3.0, 3.0};
    std::swap(states[0], states[1]);
    states[1].reset();
    EXPECT_EQ(states[0].get(), (Value{3.0, 3.0, 3.0, 3.0}));
    const auto* held   = &states[0].get();
    states[0].mutate() = Value{4.0, 4.0, 4.0, 4.0};
    EXPECT_EQ(&states[0].get(), held);

    // Save, then change the saved state's copy twice: the second write does not clone again, and the state below
    // keeps its value.
    states[1]             = states[0];
    states[1].mutate()[0] = 5.0;
    const auto* copied    = &states[1].get();
    states[1].mutate()[1] = 6.0;
    EXPECT_EQ(&states[1].get(), copied);
    EXPECT_EQ(states[0].get(), (Value{4.0, 4.0, 4.0, 4.0}));
    EXPECT_EQ(states[1].get(), (Value{5.0, 6.0, 4.0, 4.0}));
}
