//
// Dependency wiring, with a fake term and Parameter inputs and no Model: a dependent is marked dirty by a change after
// its own bookkeeping, follows its inputs' saved states, comes back from a restore with the value it had at the save,
// and settles on an accept. The ids are the Proposal's: nested saves restore innermost first.
//

#include "gtest/gtest.h"

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"
#include "core/abstract/observables/Dependency.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/computation/Computation.h"
#include "core/parameters/Parameter.h"

#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core;
using Parameter = parameters::Parameter<double>;

namespace {
    // A term summing its inputs, which records what happens to it.
    class Sum : public computation::Computation<double>,
                public abstract::Observable<Sum>,
                public abstract::Cacheable<Sum>,
                public abstract::Checkpointable<Sum, double> {
    public:
        std::vector<std::string> log;
        int evaluations = 0;

        void add(const std::shared_ptr<Parameter>& input) {
            inputs_.push_back(input);
            abstract::dependsOn(this, input, [this]() { log.emplace_back("changed"); });
        }
        void addSilently(const std::shared_ptr<Parameter>& input) {
            inputs_.push_back(input);
            abstract::savesWith(this, input);
        }

        double value() override {
            if (isDirty()) {
                ++evaluations;
                value_ = 0;
                for (const auto& input : inputs_) {
                    value_ += input->value();
                }
                setClean();
            }
            return value_;
        }

    private:
        friend class abstract::Cacheable<Sum>;
        friend class abstract::Checkpointable<Sum, double>;
        std::vector<std::shared_ptr<Parameter>> inputs_;
    };

    // A container of Parameters that passes their events through.
    class Bag : public abstract::Observable<Bag>,
                public abstract::UncacheablePassthrough<Bag>,
                public abstract::CheckpointablePassthrough<Bag> {
    public:
        explicit Bag(const std::shared_ptr<Parameter>& part) {
            abstract::passThrough(this, part);
        }
    };

    // Changes `parameter` to `value` within saved state `id`, as a Proposal does.
    void change(Parameter& parameter, const int id, const double value) {
        parameter.saveState(id);
        parameter.setValue(value);
    }
}// namespace

// A change runs the dependent's bookkeeping, then marks it dirty; it is not evaluated until read.
TEST(DependencyTest, AChangeRunsTheBookkeepingThenMarksTheDependentDirty) {
    auto x = std::make_shared<Parameter>(1.0);
    Sum sum;
    sum.add(x);
    sum.setDirty();
    EXPECT_EQ(sum.value(), 1.0);
    EXPECT_FALSE(sum.isDirty());

    change(*x, 1, 5.0);
    EXPECT_TRUE(sum.isDirty());
    EXPECT_EQ(sum.log, std::vector<std::string>{"changed"});
    EXPECT_EQ(sum.evaluations, 1);
    EXPECT_EQ(sum.value(), 5.0);
    EXPECT_EQ(sum.evaluations, 2);
}

// A restore brings the dependent back to the value it had when its input was saved, clean; an accept keeps what it
// was evaluated to.
TEST(DependencyTest, ARestoreBringsBackTheValueAtTheSaveAndAnAcceptKeepsTheNewOne) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    Sum sum;
    sum.add(x);
    sum.add(y);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 11.0);

    change(*x, 1, 2.0);
    EXPECT_EQ(sum.value(), 12.0);
    x->restoreState(1);
    EXPECT_FALSE(sum.isDirty());
    EXPECT_EQ(sum.value(), 11.0);
    EXPECT_EQ(sum.evaluations, 2);

    change(*y, 2, 20.0);
    EXPECT_EQ(sum.value(), 21.0);
    y->acceptState();
    EXPECT_FALSE(sum.isDirty());
    EXPECT_EQ(sum.value(), 21.0);
    EXPECT_EQ(sum.evaluations, 3);
}

// Saves nest by id: an inner save restores on its own, and an outer restore then undoes the rest.
TEST(DependencyTest, NestedSavesRestoreInnermostFirst) {
    auto x = std::make_shared<Parameter>(1.0);
    Sum sum;
    sum.add(x);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 1.0);

    change(*x, 1, 2.0);
    ASSERT_EQ(sum.value(), 2.0);
    change(*x, 2, 3.0);
    ASSERT_EQ(sum.value(), 3.0);

    x->restoreState(1);// not on top: ignored
    EXPECT_EQ(sum.value(), 3.0);
    x->restoreState(2);
    EXPECT_EQ(sum.value(), 2.0);
    x->restoreState(1);
    EXPECT_EQ(sum.value(), 1.0);
}

// A second save under the same id, from another input, saves nothing more: the restore is to the first save.
TEST(DependencyTest, ASecondSaveUnderTheSameIdSavesNothingMore) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    Sum sum;
    sum.add(x);
    sum.add(y);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 11.0);

    change(*x, 1, 2.0);// saves the sum at 11; the sum is dirty and unread
    change(*y, 1, 20.0);
    EXPECT_EQ(sum.evaluations, 1);
    EXPECT_EQ(sum.value(), 22.0);
    y->restoreState(1);
    x->restoreState(1);
    EXPECT_EQ(sum.value(), 11.0);
    EXPECT_EQ(sum.evaluations, 2);
}

// An accept settles the saved states and nothing else: a dependent nothing read since the change stays dirty, and
// evaluates when next read.
TEST(DependencyTest, AnAcceptLeavesAnUnreadDependentDirty) {
    auto x = std::make_shared<Parameter>(1.0);
    Sum sum;
    sum.add(x);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 1.0);

    change(*x, 1, 2.0);
    x->acceptState();
    EXPECT_TRUE(sum.isDirty());
    EXPECT_EQ(sum.value(), 2.0);
    EXPECT_FALSE(sum.isDirty());
}

// A save evaluates the dependent first, so a dependent dirty at the save restores to the value its inputs had then.
TEST(DependencyTest, ASaveEvaluatesTheDependentSoItsRestoreIsCurrent) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    Sum sum;
    sum.add(x);
    sum.add(y);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 11.0);

    change(*x, 1, 2.0);
    x->acceptState();// accepted unread: the sum is dirty
    change(*y, 2, 20.0);// y's save evaluates the sum at x = 2, y = 10 before y changes
    EXPECT_EQ(sum.evaluations, 2);
    EXPECT_EQ(sum.value(), 22.0);
    y->restoreState(2);
    EXPECT_FALSE(sum.isDirty());
    EXPECT_EQ(sum.value(), 12.0);
    EXPECT_EQ(sum.evaluations, 3);
}

// A dependent that only saves with an input is not marked dirty by its changes.
TEST(DependencyTest, SavesWithFollowsSavedStatesWithoutDirtying) {
    auto x = std::make_shared<Parameter>(1.0);
    Sum sum;
    sum.addSilently(x);
    sum.setDirty();
    ASSERT_EQ(sum.value(), 1.0);

    change(*x, 1, 2.0);
    EXPECT_FALSE(sum.isDirty());
    EXPECT_EQ(sum.value(), 1.0);
    sum.setDirty();
    EXPECT_EQ(sum.value(), 2.0);
    x->restoreState(1);
    EXPECT_EQ(sum.value(), 1.0);
}

// A container passes its part's changes and saved states through, so a dependent of the container depends on the part.
TEST(DependencyTest, PassThroughForwardsChangesAndSavedStates) {
    auto part = std::make_shared<Parameter>(1.0);
    Bag bag(part);
    std::vector<std::string> seen;
    bag.add_pre_change_listener([&]() { seen.emplace_back("pre"); });
    bag.add_post_change_listener([&]() { seen.emplace_back("post"); });
    bag.add_save_state_listener([&](const int id) { seen.push_back("save " + std::to_string(id)); });
    bag.add_accept_state_listener([&]() { seen.emplace_back("accept"); });
    bag.add_restore_state_listener([&](const int id) { seen.push_back("restore " + std::to_string(id)); });

    change(*part, 7, 2.0);
    part->restoreState(7);
    change(*part, 8, 3.0);
    part->acceptState();
    EXPECT_EQ(seen, (std::vector<std::string>{"save 7", "pre", "post", "restore 7", "save 8", "pre", "post", "accept"}));
}
