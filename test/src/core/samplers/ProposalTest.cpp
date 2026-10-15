//
// A Proposal is one Metropolis-Hastings step: it saves what it changes, evaluates the target before it settles, keeps
// or restores its changes, and counts once. Tries nest inside a move and revert without counting.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"
#include "core/samplers/Proposal.h"

#include <boost/random.hpp>

#include <limits>
#include <memory>
#include <stdexcept>

using transmission_nets::core::parameters::Parameter;
using transmission_nets::core::samplers::Proposal;
using transmission_nets::core::samplers::SamplerCounts;

namespace {
    constexpr double INF = std::numeric_limits<double>::infinity();

    // A target that is the sum of two parameters and counts how often it is evaluated.
    struct SumTarget {
        std::shared_ptr<Parameter<double>> a = std::make_shared<Parameter<double>>(1.0);
        std::shared_ptr<Parameter<double>> b = std::make_shared<Parameter<double>>(10.0);
        int evaluations                       = 0;

        double value() {
            ++evaluations;
            return a->value() + b->value();
        }
    };

    using Engine = boost::random::mt19937;

    struct Fixture {
        SumTarget target{};
        Engine rng{7};
        SamplerCounts counts{};

        Proposal<SumTarget, Engine> propose() { return Proposal<SumTarget, Engine>(target, rng, counts); }
    };
}// namespace

TEST(ProposalTest, ReadsTheCurrentValueWhenOpened) {
    Fixture f;
    auto proposal = f.propose();
    EXPECT_EQ(proposal.current(), 11.0);
    EXPECT_EQ(f.target.evaluations, 1);
    proposal.reject();
}

TEST(ProposalTest, AcceptKeepsEveryChangeAndCountsOneAccepted) {
    Fixture f;
    auto proposal = f.propose();
    proposal.change(*f.target.a, 2.0);
    proposal.change(*f.target.b, 20.0);
    EXPECT_EQ(proposal.value(), 22.0);
    proposal.accept();

    EXPECT_EQ(f.target.a->value(), 2.0);
    EXPECT_EQ(f.target.b->value(), 20.0);
    EXPECT_EQ(f.counts, (SamplerCounts{.accepted = 1, .proposed = 1}));
    // An accepted state has nothing saved to restore.
    f.target.a->restoreState(1);
    EXPECT_EQ(f.target.a->value(), 2.0);
}

TEST(ProposalTest, RejectRestoresEveryChangeAndCountsOneRejected) {
    Fixture f;
    auto proposal = f.propose();
    proposal.change(*f.target.a, 2.0);
    proposal.change(*f.target.b, 20.0);
    proposal.change(*f.target.a, 3.0);// a second change to the same parameter
    proposal.reject();

    EXPECT_EQ(f.target.a->value(), 1.0);
    EXPECT_EQ(f.target.b->value(), 10.0);
    EXPECT_EQ(f.counts, (SamplerCounts{.accepted = 0, .proposed = 1}));
}

TEST(ProposalTest, RevertRestoresWithoutCounting) {
    Fixture f;
    auto trial = f.propose();
    trial.change(*f.target.a, 2.0);
    EXPECT_EQ(trial.value(), 12.0);
    trial.revert();

    EXPECT_EQ(f.target.a->value(), 1.0);
    EXPECT_EQ(f.counts, SamplerCounts{});
}

TEST(ProposalTest, AcceptAndRejectEvaluateAnUnevaluatedChangeFirst) {
    {
        Fixture f;
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        proposal.accept();
        EXPECT_EQ(f.target.evaluations, 2);// opened, then before the accept
    }
    {
        Fixture f;
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        proposal.reject();
        EXPECT_EQ(f.target.evaluations, 2);// opened, then before the restore
    }
    {
        Fixture f;
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        proposal.value();
        proposal.reject();
        EXPECT_EQ(f.target.evaluations, 2);// already evaluated: not again, and not after the restore
    }
}

TEST(ProposalTest, AnUntouchedProposalCountsWithoutEvaluatingAgain) {
    Fixture f;
    f.propose().reject();
    f.propose().accept();
    EXPECT_EQ(f.counts, (SamplerCounts{.accepted = 1, .proposed = 2}));
    EXPECT_EQ(f.target.evaluations, 2);
}

TEST(ProposalTest, DecideAcceptsOnACertainRatioAndRejectsOnAnImpossibleOne) {
    Fixture f;
    {
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        EXPECT_TRUE(proposal.decide(INF));
        EXPECT_EQ(f.target.a->value(), 2.0);
    }
    {
        auto proposal = f.propose();
        proposal.change(*f.target.a, 3.0);
        EXPECT_FALSE(proposal.decide(-INF));
        EXPECT_EQ(f.target.a->value(), 2.0);
    }
    EXPECT_EQ(f.counts, (SamplerCounts{.accepted = 1, .proposed = 2}));
}

TEST(ProposalTest, MetropolisUsesTheRatioOfProposedToCurrentPlusTheAdjustment) {
    Fixture f;
    {
        // Proposed 12 vs current 11, adjustment makes the log ratio +inf: accepted.
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        EXPECT_TRUE(proposal.metropolis(INF));
    }
    {
        auto proposal = f.propose();
        proposal.change(*f.target.a, 3.0);
        EXPECT_FALSE(proposal.metropolis(-INF));
    }
    EXPECT_EQ(f.target.a->value(), 2.0);
}

// A multiple-try move evaluates tries while its own change is saved; each try restores only its own change.
TEST(ProposalTest, ANestedTryRestoresOnlyItsOwnChange) {
    Fixture f;
    auto move = f.propose();
    move.change(*f.target.a, 2.0);
    EXPECT_EQ(move.value(), 12.0);
    {
        auto trial = f.propose();
        trial.change(*f.target.a, 5.0);
        EXPECT_EQ(trial.value(), 15.0);
        trial.revert();
    }
    EXPECT_EQ(f.target.a->value(), 2.0);
    move.reject();
    EXPECT_EQ(f.target.a->value(), 1.0);
    EXPECT_EQ(f.counts, (SamplerCounts{.accepted = 0, .proposed = 1}));
}

TEST(ProposalTest, AChangeInsideANestedTryMakesTheMoveEvaluateAgainBeforeItSettles) {
    Fixture f;
    auto move = f.propose();
    move.change(*f.target.a, 2.0);
    move.value();
    const int evaluated = f.target.evaluations;
    {
        auto trial = f.propose();
        trial.change(*f.target.b, 50.0);
        trial.value();
        trial.revert();
    }
    const int afterTrial = f.target.evaluations;
    move.accept();
    EXPECT_EQ(f.target.evaluations, afterTrial + 1) << "the move's own evaluation predates the try's change";
    EXPECT_GT(afterTrial, evaluated);
    EXPECT_EQ(f.target.a->value(), 2.0);
    EXPECT_EQ(f.target.b->value(), 10.0);
}

TEST(ProposalTest, AnExceptionMidStepPutsEveryChangeBackWithoutCounting) {
    Fixture f;
    try {
        auto proposal = f.propose();
        proposal.change(*f.target.a, 2.0);
        proposal.change(*f.target.b, 20.0);
        throw std::runtime_error("the step failed");
    } catch (const std::runtime_error&) {
    }
    EXPECT_EQ(f.target.a->value(), 1.0);
    EXPECT_EQ(f.target.b->value(), 10.0);
    EXPECT_EQ(f.counts, SamplerCounts{});
    // The thread's open-step list is back to empty: a new step opens and settles normally.
    auto next = f.propose();
    next.change(*f.target.a, 3.0);
    next.accept();
    EXPECT_EQ(f.target.a->value(), 3.0);
}
