//
// The Accumulator's running sum: it takes a term's contribution out when the term goes dirty and adds the new one back
// when read, and its sum follows the terms' saved states, so after a change, an evaluation, a reject or an accept, and
// through a nested save, it equals the sum of the terms' values. No Model.
//

#include "gtest/gtest.h"

#include "core/abstract/observables/Dependency.h"
#include "core/computation/Accumulator.h"
#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core;
using Parameter = parameters::Parameter<double>;

namespace {
    // A term whose value is its parameter times a factor.
    class Scaled : public computation::PartialLikelihood {
    public:
        Scaled(std::shared_ptr<Parameter> x, const double factor) : x_(std::move(x)), factor_(factor) {
            abstract::dependsOn(this, x_);
            setDirty();
        }
        computation::Likelihood value() override {
            if (isDirty()) {
                ++evaluations;
                value_ = factor_ * x_->value();
                setClean();
            }
            return value_;
        }
        std::string identifier() override { return "Scaled"; }
        int evaluations = 0;

    private:
        std::shared_ptr<Parameter> x_;
        double factor_;
    };

    using Sum = computation::Accumulator<computation::PartialLikelihood, computation::Likelihood>;

    void change(Parameter& parameter, const int id, const double value) {
        parameter.saveState(id);
        parameter.setValue(value);
    }
}// namespace

// Only the terms that went dirty are evaluated again, and the sum is the sum of the terms.
TEST(AccumulatorTest, ADirtyTermIsTakenOutAndItsNewValueAddedBack) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    auto a = std::make_shared<Scaled>(x, 2.0);
    auto b = std::make_shared<Scaled>(y, 3.0);
    Sum sum;
    sum.addTarget(a);
    sum.addTarget(b);
    ASSERT_EQ(sum.value(), 32.0);

    change(*x, 1, 4.0);
    EXPECT_TRUE(sum.isDirty());
    EXPECT_EQ(sum.value(), 38.0);
    EXPECT_EQ(a->evaluations, 2);
    EXPECT_EQ(b->evaluations, 1);
    x->acceptState();
    EXPECT_EQ(sum.value(), 38.0);
}

// A rejected change restores the sum along with the term, and the restored term's value is what the next change
// takes out.
TEST(AccumulatorTest, ARejectedChangeRestoresTheSumAndTheContributionItTakesOutNext) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    auto a = std::make_shared<Scaled>(x, 2.0);
    auto b = std::make_shared<Scaled>(y, 3.0);
    Sum sum;
    sum.addTarget(a);
    sum.addTarget(b);
    ASSERT_EQ(sum.value(), 32.0);

    change(*x, 1, 4.0);
    ASSERT_EQ(sum.value(), 38.0);
    x->restoreState(1);
    EXPECT_FALSE(sum.isDirty());
    EXPECT_EQ(sum.value(), 32.0);

    change(*x, 2, 5.0);
    EXPECT_EQ(sum.value(), 40.0);
    x->acceptState();
    change(*y, 3, 1.0);
    EXPECT_EQ(sum.value(), 13.0);
    y->restoreState(3);
    EXPECT_EQ(sum.value(), 40.0);
}

// A nested change, evaluated and reverted, leaves the enclosing change's sum as it was.
TEST(AccumulatorTest, ANestedChangeRevertedLeavesTheEnclosingSum) {
    auto x = std::make_shared<Parameter>(1.0);
    auto y = std::make_shared<Parameter>(10.0);
    auto a = std::make_shared<Scaled>(x, 2.0);
    auto b = std::make_shared<Scaled>(y, 3.0);
    Sum sum;
    sum.addTarget(a);
    sum.addTarget(b);
    ASSERT_EQ(sum.value(), 32.0);

    change(*x, 1, 4.0);
    ASSERT_EQ(sum.value(), 38.0);
    change(*y, 2, 20.0);
    ASSERT_EQ(sum.value(), 68.0);
    y->restoreState(2);
    EXPECT_EQ(sum.value(), 38.0);
    x->acceptState();
    EXPECT_EQ(sum.value(), 38.0);
    change(*y, 3, 0.0);
    EXPECT_EQ(sum.value(), 8.0);
}

// A term whose value is NaN makes the sum -infinity rather than NaN.
TEST(AccumulatorTest, ANaNTermMakesTheSumMinusInfinity) {
    auto x = std::make_shared<Parameter>(std::numeric_limits<double>::quiet_NaN());
    auto a = std::make_shared<Scaled>(x, 1.0);
    Sum sum;
    sum.addTarget(a);
    EXPECT_EQ(sum.value(), -std::numeric_limits<double>::infinity());
}
