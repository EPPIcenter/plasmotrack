//
// ReplicaExchange reports its ensemble through queries and return values: what occupies each temperature slot, how
// often neighbouring slots swapped, and what temperature adaptation decided. It records everything it carries between
// Sampling iterations and restores from that record; a Checkpoint is that record on disk.
//

#include "gtest/gtest.h"

#include "core/samplers/meta/ReplicaExchange.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/EnsembleCheckpoint.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/SamplerSet.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::samplers::ReplicaExchange;

namespace {
    using Exchange = ReplicaExchange<State, Model, SamplerSet>;

    nlohmann::json exampleInput() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return nlohmann::json::parse(in);
    }

    std::unique_ptr<Exchange> makeExchange(const int numChains, const unsigned int seed) {
        const std::vector<double> idp(60, 1.0 / 60);
        auto rng = std::make_shared<boost::random::mt19937>(seed);
        return std::make_unique<Exchange>(numChains, 0.5, rng, false, 1, exampleInput(), idp, idp);
    }
}// namespace

TEST(ReplicaExchangeTest, SlotsStartInReplicaOrderOnAnEvenlySpacedLadder) {
    auto repex = makeExchange(4, 1);

    const auto slots = repex->slots();

    ASSERT_EQ(slots.size(), 4u);
    const std::vector<double> expected{1.0, 0.875, 0.75, 0.625};
    for (size_t slot = 0; slot < slots.size(); ++slot) {
        EXPECT_EQ(slots[slot].replica, slot);
        EXPECT_DOUBLE_EQ(slots[slot].temperature, expected[slot]);
    }
}

TEST(ReplicaExchangeTest, SlotZeroIsTheHotChainWithItsPriorAndLikelihood) {
    auto repex = makeExchange(4, 2);
    for (int ii = 0; ii < 30; ++ii) {
        repex->burnin();
    }

    const auto slots = repex->slots();
    const auto& hot  = repex->hotChain();

    EXPECT_EQ(&repex->chains[slots[0].replica], &hot);
    EXPECT_DOUBLE_EQ(slots[0].temperature, 1.0);
    EXPECT_DOUBLE_EQ(slots[0].prior, hot.model->getPrior());
    EXPECT_DOUBLE_EQ(slots[0].likelihood, hot.model->getLikelihood());
    for (const auto& slot : slots) {
        const auto& chain = repex->chains[slot.replica];
        EXPECT_DOUBLE_EQ(slot.temperature, chain.model->getTemperature());
        EXPECT_DOUBLE_EQ(slot.likelihood, chain.model->getLikelihood());
    }
}

TEST(ReplicaExchangeTest, SwapsAreAttemptedOnAlternatingPairsInBothPhases) {
    auto repex = makeExchange(4, 3);
    for (int ii = 0; ii < 3; ++ii) {
        repex->burnin();
    }
    for (int ii = 0; ii < 2; ++ii) {
        repex->sample();
    }

    // Iterations alternate between pairs (0,1),(2,3) and pair (1,2), starting with the first.
    const auto counts = repex->swapCounts();
    ASSERT_EQ(counts.size(), 3u);
    EXPECT_EQ(counts[0].attempted, 3u);
    EXPECT_EQ(counts[1].attempted, 2u);
    EXPECT_EQ(counts[2].attempted, 3u);
}

TEST(ReplicaExchangeTest, AcceptedSwapsCountEveryChangeOfHotChain) {
    auto repex = makeExchange(2, 4);
    unsigned int handoffs = 0;
    size_t hot = repex->slots()[0].replica;
    for (int ii = 0; ii < 200; ++ii) {
        if (ii < 100) {
            repex->burnin();
        } else {
            repex->sample();
        }
        const size_t now = repex->slots()[0].replica;
        handoffs += now != hot;
        hot = now;
    }

    const auto counts = repex->swapCounts();
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_GT(handoffs, 0u);
    EXPECT_EQ(counts[0].accepted, handoffs);
    EXPECT_EQ(counts[0].attempted, 100u);
}

TEST(ReplicaExchangeTest, AdaptationReturnsTheLadderItApplied) {
    auto repex = makeExchange(4, 5);
    for (int ii = 0; ii < 15; ++ii) {
        ASSERT_FALSE(repex->burnin());
    }

    const auto adaptation = repex->burnin();

    ASSERT_TRUE(adaptation);
    ASSERT_TRUE(adaptation->applied);
    ASSERT_EQ(adaptation->temperatures.size(), 4u);
    EXPECT_DOUBLE_EQ(adaptation->temperatures.front(), 1.0);
    EXPECT_DOUBLE_EQ(adaptation->temperatures.back(), 0.625);
    for (size_t slot = 1; slot < adaptation->temperatures.size(); ++slot) {
        EXPECT_LT(adaptation->temperatures[slot], adaptation->temperatures[slot - 1]);
    }
    const auto slots = repex->slots();
    for (size_t slot = 0; slot < slots.size(); ++slot) {
        EXPECT_DOUBLE_EQ(slots[slot].temperature, adaptation->temperatures[slot]);
    }
    ASSERT_EQ(adaptation->barriers.size(), 3u);
    for (const double barrier : adaptation->barriers) {
        EXPECT_GE(barrier, 0.0);
        EXPECT_LE(barrier, 1.0);
    }
}

namespace {
    // Burn in past the last temperature adaptation, so swap barriers are pending, then sample.
    std::unique_ptr<Exchange> runPastAnAdaptation(const unsigned int seed) {
        auto repex = makeExchange(4, seed);
        for (int ii = 0; ii < 43; ++ii) {
            repex->burnin();
        }
        for (int ii = 0; ii < 15; ++ii) {
            repex->sample();
        }
        return repex;
    }
}// namespace

TEST(ReplicaExchangeTest, ExchangeStateRestoresSlotsLadderAndCounts) {
    const auto source   = runPastAnAdaptation(7);
    const auto saved    = source->exchangeState();
    auto restored       = makeExchange(4, 8);

    restored->restoreExchangeState(saved);

    EXPECT_EQ(restored->exchangeState(), saved);
    const auto expected = source->slots();
    const auto actual   = restored->slots();
    for (size_t slot = 0; slot < expected.size(); ++slot) {
        EXPECT_EQ(actual[slot].replica, expected[slot].replica);
        EXPECT_DOUBLE_EQ(actual[slot].temperature, expected[slot].temperature);
    }
    EXPECT_EQ(&restored->hotChain(), &restored->chains[expected[0].replica]);
    for (size_t pair = 0; pair < 3; ++pair) {
        EXPECT_EQ(restored->swapCounts()[pair].accepted, source->swapCounts()[pair].accepted);
        EXPECT_EQ(restored->swapCounts()[pair].attempted, source->swapCounts()[pair].attempted);
    }
}

TEST(ReplicaExchangeTest, RestoredExchangeAlternatesAsTheSource) {
    const auto source = runPastAnAdaptation(9);
    auto restored     = makeExchange(4, 10);
    restored->restoreExchangeState(source->exchangeState());

    // The next iteration attempts the same pairs in both.
    source->sample();
    restored->sample();
    for (size_t pair = 0; pair < 3; ++pair) {
        EXPECT_EQ(restored->swapCounts()[pair].attempted, source->swapCounts()[pair].attempted);
    }
}

TEST(ReplicaExchangeTest, RestoreRefusesAStateForAnotherEnsemble) {
    const auto saved = runPastAnAdaptation(11)->exchangeState();
    auto twoChains   = makeExchange(2, 12);
    EXPECT_THROW(twoChains->restoreExchangeState(saved), std::invalid_argument);

    auto repeated     = saved;
    repeated.replicas = {0, 0, 1, 2};
    auto fourChains   = makeExchange(4, 12);
    EXPECT_THROW(fourChains->restoreExchangeState(repeated), std::invalid_argument);

    auto shortLadder = saved;
    shortLadder.temperatures.pop_back();
    EXPECT_THROW(fourChains->restoreExchangeState(shortLadder), std::invalid_argument);
}

namespace {
    // The burn-in iterations, counted from 0, among `from` to `to` (exclusive) whose burnin() call adapted the ladder.
    template<typename AnyExchange>
    std::vector<int> adaptingIterations(AnyExchange& repex, const int from, const int to) {
        std::vector<int> adapted;
        for (int ii = from; ii < to; ++ii) {
            if (repex.burnin()) {
                adapted.push_back(ii);
            }
        }
        return adapted;
    }
}// namespace

TEST(ReplicaExchangeTest, LadderAdaptsAtBurnInIterationFifteenThenEveryFive) {
    auto repex = makeExchange(4, 13);

    EXPECT_EQ(adaptingIterations(*repex, 0, 31), (std::vector<int>{15, 20, 25, 30}));
}

TEST(ReplicaExchangeTest, LadderNeverAdaptsWithOneOrTwoChains) {
    for (const int numChains : {1, 2}) {
        auto repex = makeExchange(numChains, 14);

        EXPECT_EQ(adaptingIterations(*repex, 0, 31), std::vector<int>{}) << numChains << " chains";
    }
}

TEST(ReplicaExchangeTest, RestoredExchangeAdaptsAtTheSameBurnInIterations) {
    auto source = makeExchange(4, 15);
    adaptingIterations(*source, 0, 17);
    auto restored = makeExchange(4, 16);

    restored->restoreExchangeState(source->exchangeState());

    EXPECT_EQ(adaptingIterations(*restored, 17, 31), (std::vector<int>{20, 25, 30}));
}

TEST(ReplicaExchangeTest, SnapshotMatchesTheSlotSwapAndSamplerCountQueries) {
    auto repex = makeExchange(4, 17);
    for (int ii = 0; ii < 20; ++ii) {
        repex->burnin();
    }
    for (int ii = 0; ii < 5; ++ii) {
        repex->sample();
    }

    const auto snapshot = repex->snapshot();

    EXPECT_EQ(&snapshot.hot.state, repex->hotChain().state.get());
    EXPECT_EQ(&snapshot.hot.model, repex->hotChain().model.get());
    const auto slots = repex->slots();
    ASSERT_EQ(snapshot.slots.size(), slots.size());
    for (size_t slot = 0; slot < slots.size(); ++slot) {
        EXPECT_EQ(snapshot.slots[slot].replica, slots[slot].replica);
        EXPECT_EQ(snapshot.slots[slot].temperature, slots[slot].temperature);
        EXPECT_EQ(snapshot.slots[slot].prior, slots[slot].prior);
        EXPECT_EQ(snapshot.slots[slot].likelihood, slots[slot].likelihood);
    }
    EXPECT_EQ(snapshot.swaps, repex->swapCounts());

    const auto& ids = repex->samplerIds();
    EXPECT_EQ(ids, repex->chains.front().sampler->ids());
    ASSERT_EQ(snapshot.acceptance.size(), 4u);
    for (size_t replica = 0; replica < 4; ++replica) {
        ASSERT_EQ(snapshot.acceptance[replica].size(), ids.size());
        for (size_t id = 0; id < ids.size(); ++id) {
            EXPECT_EQ(snapshot.acceptance[replica][id], repex->chains[replica].sampler->sampler(ids[id]).counts()) << ids[id];
        }
    }
}

TEST(ReplicaExchangeTest, FinalAdaptationIsTheHotChainsSamplersAndTheLadder) {
    auto repex = makeExchange(4, 18);
    for (int ii = 0; ii < 30; ++ii) {
        repex->burnin();
    }

    const auto adaptation = repex->finalAdaptation();

    const auto& ids       = repex->samplerIds();
    const auto& samplers  = *repex->hotChain().sampler;
    ASSERT_EQ(adaptation.samplers.size(), ids.size());
    for (size_t id = 0; id < ids.size(); ++id) {
        EXPECT_EQ(adaptation.samplers[id], samplers.sampler(ids[id]).adaptationState()) << ids[id];
    }
    const auto slots = repex->slots();
    ASSERT_EQ(adaptation.temperatures.size(), slots.size());
    for (size_t slot = 0; slot < slots.size(); ++slot) {
        EXPECT_EQ(adaptation.temperatures[slot], slots[slot].temperature);
    }
}

namespace {
    // A replica that never moves: its log-likelihood is fixed, so the swap acceptance between two slots is known.
    struct FixedState {
        double logLikelihood;
    };

    struct FixedModel {
        FixedModel(std::shared_ptr<FixedState> state, const double temperature) : state(std::move(state)), temperature(temperature) {}
        double getLikelihood() const { return state->logLikelihood; }
        double getPrior() const { return 0.0; }
        double getTemperature() const { return temperature; }
        void setTemperature(const double t) { temperature = t; }

        std::shared_ptr<FixedState> state;
        double temperature;
    };

    // Registers no samplers and steps nothing.
    template<typename... Model>
    struct NoSamplers {
        NoSamplers(std::shared_ptr<FixedState>, std::shared_ptr<FixedModel>, std::shared_ptr<boost::random::mt19937>) {}
        void step(bool) {}

        const std::vector<std::string>& ids() const { return ids_; }
        transmission_nets::core::samplers::AbstractSampler& sampler(const std::string&) const { throw std::logic_error("no samplers"); }
        std::vector<std::string> ids_;
    };

    using FixedExchange = ReplicaExchange<FixedState, FixedModel, NoSamplers>;

    // Replica `r` has log-likelihood -1e6 r, so each colder slot holds a far worse replica and every swap is rejected
    // with probability 1 to double precision.
    std::unique_ptr<FixedExchange> makeRejectingExchange(const int numChains) {
        auto rng = std::make_shared<boost::random::mt19937>(20);
        return std::make_unique<FixedExchange>(numChains, 0.5, rng, 1, [](const size_t replica, auto) {
            return std::make_shared<FixedState>(FixedState{-1e6 * static_cast<double>(replica)});
        });
    }

    unsigned int acceptedSwaps(const FixedExchange& repex) {
        unsigned int accepted = 0;
        for (const auto& counts : repex.swapCounts()) {
            accepted += counts.accepted;
        }
        return accepted;
    }
}// namespace

TEST(ReplicaExchangeTest, AdaptationLeavesAnEvenlySpacedLadderWhereItIsWhenEverySwapIsRejected) {
    for (const int numChains : {4, 5, 6}) {
        auto repex          = makeRejectingExchange(numChains);
        const auto expected = repex->slots();

        for (int ii = 0; ii < 41; ++ii) {
            const auto adaptation = repex->burnin();
            ASSERT_EQ(acceptedSwaps(*repex), 0u) << numChains << " chains, burn-in iteration " << ii;
            if (!adaptation) {
                continue;
            }
            ASSERT_EQ(adaptation->temperatures.size(), expected.size());
            for (size_t slot = 0; slot < expected.size(); ++slot) {
                EXPECT_NEAR(adaptation->temperatures[slot], expected[slot].temperature, 1e-12)
                        << numChains << " chains, burn-in iteration " << ii << ", slot " << slot;
            }
        }
    }
}

TEST(ReplicaExchangeTest, BarriersAreTheMeanRejectionPerAttemptOfEachPair) {
    auto repex = makeRejectingExchange(4);
    adaptingIterations(*repex, 0, 16);
    // Four iterations since the last adaptation, two per set of pairs, with the pairs from slot 0 tried next.
    auto pending              = repex->exchangeState();
    pending.swapBarriers      = {0.6, 0.4, 1.2};
    pending.barrierIterations = 4;
    pending.burnInIterations  = 20;
    pending.secondPairSetNext = false;
    repex->restoreExchangeState(pending);

    // Iteration 20 tries pairs 0 and 2 once more, rejecting both.
    const auto adaptation = repex->burnin();

    ASSERT_TRUE(adaptation);
    ASSERT_EQ(acceptedSwaps(*repex), 0u);
    ASSERT_EQ(adaptation->barriers.size(), 3u);
    EXPECT_DOUBLE_EQ(adaptation->barriers[0], 1.6 / 3);
    EXPECT_DOUBLE_EQ(adaptation->barriers[1], 0.4 / 2);
    EXPECT_DOUBLE_EQ(adaptation->barriers[2], 2.2 / 3);
}

TEST(ReplicaExchangeTest, APairWithNoSwapAttemptsLeavesTheLadderWhereItIs) {
    auto repex = makeRejectingExchange(4);
    const auto ladder         = repex->exchangeState().temperatures;
    auto pending              = repex->exchangeState();
    pending.burnInIterations  = 15;
    repex->restoreExchangeState(pending);

    // Iteration 15 is the only one since the last adaptation, and it tries pairs 0 and 2 alone.
    const auto adaptation = repex->burnin();

    ASSERT_TRUE(adaptation);
    EXPECT_FALSE(adaptation->applied);
    EXPECT_EQ(adaptation->temperatures, ladder);
    EXPECT_EQ(repex->exchangeState().temperatures, ladder);
}

// The ensemble record: what replica exchange carries between Sampling iterations, restored into a new ensemble.
namespace {
    using ChainState = checkpoint::ChainState;
    using Record     = transmission_nets::core::samplers::EnsembleRecord<ChainState>;

    Record recordOf(const Exchange& repex) {
        return repex.record([](const State& state) { return checkpoint::captureState(state); });
    }

    // The ensemble `saved` records, its States rebuilt from the example input.
    std::unique_ptr<Exchange> restoredFrom(const Record& saved, const unsigned int seed = 99) {
        const std::vector<double> idp(60, 1.0 / 60);
        auto rng = std::make_shared<boost::random::mt19937>(seed);
        return std::make_unique<Exchange>(saved, rng, 1, [&](const size_t, const ChainState& state, std::shared_ptr<boost::random::mt19937> chainRng) {
            auto fresh = std::make_shared<State>(exampleInput(), idp, idp, std::move(chainRng), false);
            checkpoint::restoreState(checkpoint::stateLayout(*fresh), state, *fresh);
            return fresh;
        });
    }

    // Every proposal scale of every sampler of chain 0.
    std::vector<double> proposalScales(const Exchange& repex) {
        const auto& samplers = *repex.chains[0].sampler;
        std::vector<double> variances;
        for (const auto& id : samplers.ids()) {
            for (const auto& scale : samplers.sampler(id).adaptationState().scales) {
                variances.push_back(scale.variance);
            }
        }
        return variances;
    }
}// namespace

// A restored ensemble records the same: every engine continues where it was, each sampler keeps its counts and
// adaptation state, each temperature slot its replica and temperature, and each chain its State.
TEST(ReplicaExchangeTest, ARecordRestoresEnginesSamplersSlotsAndStates) {
    auto source = makeExchange(3, 23);
    adaptingIterations(*source, 0, 40);
    const auto saved = recordOf(*source);

    auto restored = restoredFrom(saved);

    EXPECT_EQ(recordOf(*restored), saved);
    EXPECT_EQ(restored->exchangeState(), source->exchangeState());
    for (size_t replica = 0; replica < 3; ++replica) {
        for (int draw = 0; draw < 5; ++draw) {
            EXPECT_EQ((*restored->chains[replica].r)(), (*source->chains[replica].r)()) << "replica " << replica;
        }
        EXPECT_NEAR(restored->chains[replica].model->getLikelihood(), source->chains[replica].model->getLikelihood(), 1e-8 * std::abs(source->chains[replica].model->getLikelihood()));
    }
}

TEST(ReplicaExchangeTest, ARecordOfAnotherEnsembleIsRefused) {
    const auto saved = recordOf(*makeExchange(2, 24));

    auto otherSamplers = saved;
    otherSamplers.samplerIds.pop_back();
    EXPECT_THROW(restoredFrom(otherSamplers), std::invalid_argument);

    auto fewerCounts = saved;
    fewerCounts.replicas[0].samplers.counts.pop_back();
    EXPECT_THROW(restoredFrom(fewerCounts), std::invalid_argument);

    auto noReplicas = saved;
    noReplicas.replicas.clear();
    EXPECT_THROW(restoredFrom(noReplicas), std::invalid_argument);
}

// Samplers adapt from Sampling iteration 20 to 200 of their chain; restored at 199 they adapt once more, then stop.
TEST(ReplicaExchangeTest, RestoredSamplersAdaptOnTheSameSamplingIterations) {
    auto source = makeExchange(1, 15);
    adaptingIterations(*source, 0, 199);
    auto restored     = restoredFrom(recordOf(*source));
    const auto atSave = proposalScales(*restored);

    restored->burnin();
    const auto afterIteration199 = proposalScales(*restored);
    restored->burnin();

    EXPECT_NE(afterIteration199, atSave);
    EXPECT_EQ(proposalScales(*restored), afterIteration199);
}

TEST(ReplicaExchangeTest, ARestoredEnsembleAdaptsItsLadderIdentically) {
    // Three iterations past the adaptation at 15, so the next one's window straddles the record, with an odd count.
    auto source = makeExchange(4, 21);
    adaptingIterations(*source, 0, 18);
    auto restored = restoredFrom(recordOf(*source));

    for (int ii = 18; ii < 31; ++ii) {
        const auto expected = source->burnin();
        const auto actual   = restored->burnin();
        ASSERT_EQ(actual.has_value(), expected.has_value()) << "burn-in iteration " << ii;
        if (!expected) {
            continue;
        }
        // Restore recomputes each likelihood from its State, which can differ from the running sum in the last bits.
        EXPECT_EQ(actual->applied, expected->applied) << "burn-in iteration " << ii;
        ASSERT_EQ(actual->temperatures.size(), expected->temperatures.size());
        for (size_t slot = 0; slot < expected->temperatures.size(); ++slot) {
            EXPECT_NEAR(actual->temperatures[slot], expected->temperatures[slot], 1e-12) << "burn-in iteration " << ii << ", slot " << slot;
        }
        ASSERT_EQ(actual->barriers.size(), expected->barriers.size());
        for (size_t pair = 0; pair < expected->barriers.size(); ++pair) {
            EXPECT_NEAR(actual->barriers[pair], expected->barriers[pair], 1e-12) << "burn-in iteration " << ii << ", pair " << pair;
        }
    }
}

// A restarted burn-in keeps the State, the ladder, the slots and the swap counts, counts burn-in from 0 and has no
// pending barriers, so the ladder adapts on the iterations a fresh Run's does.
TEST(ReplicaExchangeTest, ARestartedBurnInCountsFromZeroWithNoPendingBarriers) {
    auto repex = makeExchange(4, 22);
    adaptingIterations(*repex, 0, 18);
    const auto before = repex->exchangeState();
    ASSERT_GT(before.barrierIterations, 0);

    repex->restartBurnIn();

    const auto after = repex->exchangeState();
    EXPECT_EQ(after.burnInIterations, 0);
    EXPECT_EQ(after.barrierIterations, 0);
    EXPECT_EQ(after.swapBarriers, (std::vector<double>{0.0, 0.0, 0.0}));
    EXPECT_EQ(after.temperatures, before.temperatures);
    EXPECT_EQ(after.replicas, before.replicas);
    EXPECT_EQ(after.swapCounts, before.swapCounts);
    EXPECT_EQ(adaptingIterations(*repex, 0, 31), (std::vector<int>{15, 20, 25, 30}));
}
