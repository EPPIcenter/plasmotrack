//
// The Ordering and the Candidate parents it owns: moving an Infection's duration moves it past its neighbours, each
// crossing tells only the Infections whose Candidate parents changed, the view answers the new order, and a restore
// brings it back. Four Infections observed ten days apart, every one permitted every other unless a test says
// otherwise.
//

#include "gtest/gtest.h"

#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/datatypes/Alleles.h"

#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::containers;
using transmission_nets::core::datatypes::AllelesBitSet;

namespace {
    using Genetics  = AllelesBitSet<32>;
    using Infection = transmission_nets::core::containers::Infection<Genetics>;
    using Ordering  = ObservationTimeDerivedOrdering<Infection>;
    using Ptr       = std::shared_ptr<Infection>;

    struct Fixture {
        std::vector<Ptr> infections;
        std::shared_ptr<Ordering> ordering = std::make_shared<Ordering>();
        std::vector<std::string> toldOfCandidates;// ids, in the order they were told
        std::vector<std::string> toldOfGenotypes;

        // `permits(child id, parent id)` says what the study permits; every pair by default.
        explicit Fixture(const std::function<bool(const std::string&, const std::string&)>& permits = [](const std::string&, const std::string&) { return true; }) {
            for (const auto* id : {"1", "2", "3", "4"}) {
                infections.push_back(std::make_shared<Infection>(id, 1000 + 10 * (id[0] - '1'), false));
            }
            ordering->addElements(infections, [&](const Ptr& child) {
                std::vector<Ptr> permitted;
                for (const auto& parent : infections) {
                    if (parent != child and permits(child->id(), parent->id())) {
                        permitted.push_back(parent);
                    }
                }
                return permitted;
            });
            for (const auto& infection : infections) {
                ordering->add_keyed_candidates_changed_listener(infection, [this, id = infection->id()]() { toldOfCandidates.push_back(id); });
                ordering->add_keyed_candidate_genotype_changed_listener(infection, [this, id = infection->id()]() { toldOfGenotypes.push_back(id); });
            }
        }

        Ptr operator[](const std::string& id) const {
            for (const auto& infection : infections) {
                if (infection->id() == id) {
                    return infection;
                }
            }
            throw std::out_of_range(id);
        }

        // Ids of `id`'s Candidate parents, in the view's order.
        std::vector<std::string> candidates(const std::string& id) {
            std::vector<std::string> ids;
            for (const auto& parent : ordering->candidateParents((*this)[id])) {
                ids.push_back(parent->id());
            }
            return ids;
        }

        std::vector<std::string> order() {
            std::vector<std::string> ids;
            for (const auto& infection : ordering->value()) {
                ids.push_back(infection->id());
            }
            return ids;
        }

        // Moves `id`'s infection time by `days`, under saved state 1, as a duration proposal does.
        void move(const std::string& id, const double days) {
            auto duration = (*this)[id]->infectionDuration();
            duration->saveState(1);
            duration->setValue(duration->value() - days);
        }
        void undo(const std::string& id) { (*this)[id]->infectionDuration()->restoreState(1); }
    };

    using Ids = std::vector<std::string>;
}// namespace

TEST(OrderingTest, CandidateParentsAreThePermittedInfectionsThatPrecede) {
    Fixture f;
    EXPECT_EQ(f.order(), (Ids{"1", "2", "3", "4"}));
    EXPECT_EQ(f.candidates("1"), Ids{});
    EXPECT_EQ(f.candidates("2"), (Ids{"1"}));
    EXPECT_EQ(f.candidates("3"), (Ids{"1", "2"}));
    EXPECT_EQ(f.candidates("4"), (Ids{"1", "2", "3"}));
    EXPECT_EQ(f.ordering->position(f["3"]), 2u);

    Fixture restricted([](const std::string& child, const std::string& parent) { return child != "4" or parent == "2"; });
    EXPECT_EQ(restricted.candidates("4"), (Ids{"2"}));
    EXPECT_EQ(restricted.candidates("3"), (Ids{"1", "2"}));
}

TEST(OrderingTest, AMoveWithinTheGapsChangesNothing) {
    Fixture f;
    const auto version = f.ordering->version();

    f.move("2", 5);

    EXPECT_EQ(f.order(), (Ids{"1", "2", "3", "4"}));
    EXPECT_EQ(f.ordering->version(), version);
    EXPECT_TRUE(f.toldOfCandidates.empty());
    EXPECT_EQ(f.candidates("3"), (Ids{"1", "2"}));
}

// Moving 1 past 2: 1 gains 2 as a Candidate parent and 2 loses 1; nobody else is told. A restore brings the view back
// and tells nobody, since what depends on it restores itself.
TEST(OrderingTest, AMoveTellsOnlyTheInfectionsWhoseCandidateParentsChanged) {
    Fixture f;
    const auto version = f.ordering->version();

    f.move("1", 15);

    EXPECT_EQ(f.order(), (Ids{"2", "1", "3", "4"}));
    EXPECT_EQ(f.toldOfCandidates, (Ids{"1", "2"}));
    EXPECT_GT(f.ordering->version(), version);
    EXPECT_EQ(f.candidates("1"), (Ids{"2"}));
    EXPECT_EQ(f.candidates("2"), Ids{});
    EXPECT_EQ(f.candidates("3"), (Ids{"1", "2"}));
    const auto moved = f.ordering->version();

    f.undo("1");

    EXPECT_EQ(f.order(), (Ids{"1", "2", "3", "4"}));
    EXPECT_EQ(f.toldOfCandidates, (Ids{"1", "2"}));
    EXPECT_GT(f.ordering->version(), moved);
    EXPECT_EQ(f.candidates("1"), Ids{});
    EXPECT_EQ(f.candidates("2"), (Ids{"1"}));
}

TEST(OrderingTest, AMovePastSeveralTellsEachCrossing) {
    Fixture f;

    f.move("1", 25);

    EXPECT_EQ(f.order(), (Ids{"2", "3", "1", "4"}));
    EXPECT_EQ(f.toldOfCandidates, (Ids{"1", "2", "1", "3"}));
    EXPECT_EQ(f.candidates("1"), (Ids{"2", "3"}));
    EXPECT_EQ(f.candidates("3"), (Ids{"2"}));
    EXPECT_EQ(f.candidates("4"), (Ids{"1", "2", "3"}));

    f.undo("1");
    EXPECT_EQ(f.candidates("3"), (Ids{"1", "2"}));
}

TEST(OrderingTest, AMoveEarlierTellsTheInfectionsPassed) {
    Fixture f;

    f.move("4", -25);

    EXPECT_EQ(f.order(), (Ids{"1", "4", "2", "3"}));
    EXPECT_EQ(f.toldOfCandidates, (Ids{"3", "4", "2", "4"}));
    EXPECT_EQ(f.candidates("4"), (Ids{"1"}));
    EXPECT_EQ(f.candidates("2"), (Ids{"1", "4"}));
}

// A crossing is told only to Infections the study permits the crossing pair: 4 may have no parent but 2.
TEST(OrderingTest, ACrossingTheStudyDoesNotPermitTellsNobody) {
    Fixture f([](const std::string& child, const std::string& parent) { return child != "4" or parent == "2"; });

    f.move("3", 15);

    EXPECT_EQ(f.order(), (Ids{"1", "2", "4", "3"}));
    EXPECT_EQ(f.toldOfCandidates, (Ids{"3"}));// 3 gained 4; 4 may not have 3
    EXPECT_EQ(f.candidates("4"), (Ids{"2"}));
    EXPECT_EQ(f.candidates("3"), (Ids{"1", "2", "4"}));
}

TEST(OrderingTest, RestoreOrderPutsTheElementsWhereSaidWithoutTellingAnyone) {
    Fixture f;
    const auto version = f.ordering->version();

    f.ordering->restoreOrder({f["4"], f["3"], f["2"], f["1"]});

    EXPECT_EQ(f.order(), (Ids{"4", "3", "2", "1"}));
    EXPECT_GT(f.ordering->version(), version);
    EXPECT_TRUE(f.toldOfCandidates.empty());
    EXPECT_EQ(f.candidates("1"), (Ids{"2", "3", "4"}));
    EXPECT_EQ(f.candidates("4"), Ids{});
    EXPECT_EQ(f.ordering->position(f["4"]), 0u);
}

// A genotype change of 2 reaches the Infections 2 is a Candidate parent of, in input order, and nobody before it.
TEST(OrderingTest, AGenotypeChangeIsRelayedToTheInfectionsItIsACandidateParentOf) {
    Fixture f;
    auto locus = std::make_shared<Locus>("L0", 2);
    for (const auto& infection : f.infections) {
        infection->addGenetics(locus, Genetics("10"), Genetics("10"));
    }

    auto genotype = f["2"]->latentGenotype(locus);
    genotype->saveState(1);
    genotype->setValue(Genetics("01"));

    EXPECT_EQ(f.toldOfGenotypes, (Ids{"3", "4"}));
    EXPECT_TRUE(f.toldOfCandidates.empty());
}
