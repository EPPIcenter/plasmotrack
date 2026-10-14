//
// After real MCMC sweeps (proposals, rejections, restores), every node's cached transmission likelihood
// must equal a likelihood computed from scratch. A stale or colliding parent-set cache entry breaks this.
//

#include "gtest/gtest.h"

#include "impl/model/Model/config.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/SamplerSet.h"
#include "core/samplers/Proposal.h"
#include "impl/model/Model/State.h"
#include "../model/transmission_process/node_transmission_process/MarginalLatentParentReference.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;

namespace {
    // Clonal clusters, polyclonal and missing genotypes, and allowed parents observed later than the
    // child, as in real inputs.
    nlohmann::json makeInput(boost::random::mt19937& rng) {
        const std::vector<int> alleleCounts{2, 3, 4, 5, 6, 8};
        constexpr int numNodes    = 14;
        constexpr int numClusters = 3;

        boost::random::uniform_01<> unif;
        auto randomGenotype = [&](const int numAlleles) {
            std::string bits(numAlleles, '0');
            boost::random::uniform_int_distribution<> pick(0, numAlleles - 1);
            bits[pick(rng)] = '1';
            if (unif(rng) < 0.2) {
                bits[pick(rng)] = '1';
            }
            return bits;
        };

        std::vector<std::vector<std::string>> founders(numClusters);
        for (auto& founder : founders) {
            for (const auto count : alleleCounts) {
                founder.push_back(randomGenotype(count));
            }
        }

        nlohmann::json input;
        for (size_t l = 0; l < alleleCounts.size(); ++l) {
            input["loci"].push_back({{"locus", "L" + std::to_string(l)},
                                     {"num_alleles", alleleCounts[l]},
                                     {"allele_freqs", std::vector<double>(alleleCounts[l], 1.0 / alleleCounts[l])}});
        }
        for (int n = 0; n < numNodes; ++n) {
            nlohmann::json node;
            node["id"]               = "n" + std::to_string(n);
            node["observation_time"] = 400 + 15 * n;
            node["symptomatic"]      = n % 4 != 0;
            for (size_t l = 0; l < alleleCounts.size(); ++l) {
                std::string genotype = unif(rng) < 0.15 ? randomGenotype(alleleCounts[l]) : founders[n % numClusters][l];
                if (unif(rng) < 0.05) {
                    genotype = std::string(alleleCounts[l], '0');// missing
                }
                node["observed_genotype"].push_back({{"locus", "L" + std::to_string(l)}, {"genotype", genotype}});
            }
            node["allowed_parents"] = nlohmann::json::array();
            for (int p = 0; p < numNodes; ++p) {
                if (p != n) {
                    node["allowed_parents"].push_back("n" + std::to_string(p));
                }
            }
            input["nodes"].push_back(node);
        }
        return input;
    }
}// namespace

namespace {
    void checkCacheConsistency(const std::shared_ptr<boost::random::mt19937>& rng, const nlohmann::json& input);
}// namespace

TEST(ModelCacheConsistencyTest, CachedTransmissionLikelihoodsMatchRecomputation) {
    auto rng = std::make_shared<boost::random::mt19937>(20261005);
    checkCacheConsistency(rng, makeInput(*rng));
}

namespace {
void checkCacheConsistency(const std::shared_ptr<boost::random::mt19937>& rng, const nlohmann::json& input) {
    const std::vector<double> idp(60, 1.0 / 60);


    auto state = std::make_shared<State>(input, idp, idp, rng);
    auto model = std::make_shared<Model>(state);
    SamplerSet<Model> samplers(state, model, rng);

    constexpr int sweeps = 6;
    for (int s = 0; s < sweeps; ++s) {
        samplers.step();
    }

    // MCMC acceptance is random, so also force accepted and rejected changes through every invalidation
    // path: genotypes, allele frequencies, the COI and strain-count parameters, and infection durations (which
    // reorder parent sets).
    boost::random::uniform_01<> unif;
    boost::random::uniform_int_distribution<> pickNode(0, static_cast<int>(state->infections.size()) - 1);
    auto randomBits = [&](const unsigned int numAlleles) {
        std::string bits(numAlleles, '0');
        boost::random::uniform_int_distribution<> pick(0, static_cast<int>(numAlleles) - 1);
        bits[pick(*rng)] = '1';
        if (unif(*rng) < 0.3) {
            bits[pick(*rng)] = '1';
        }
        return bits;
    };
    // Each perturbation is a Proposal, as a sampler's would be. Like a sampler, never accept an impossible state: the
    // model's running sums cannot recover from -inf.
    transmission_nets::core::samplers::SamplerCounts perturbations{};
    auto perturb = [&](auto& param, const auto& proposed) {
        auto step = transmission_nets::core::samplers::Proposal(*model, *rng, perturbations);
        step.change(*param, proposed);
        if (std::isfinite(step.value()) && unif(*rng) < 0.5) {
            step.accept();
        } else {
            step.reject();
        }
        model->value();
    };
    // Independent of the parent-set cache: evaluates every parent set directly, in the order
    // OrderBasedTransmissionProcessV3::calcParentSetDist() reports them. The process groups parent sets by their
    // parents' genomes, so its log-sum-exp agrees up to rounding.
    using ParentSetLliks = std::vector<std::pair<double, transmission_nets::core::containers::ParentSet<InfectionEvent>>>;
    // A parent set with the latent parent: its expectation over the latent parent drawn from the Source, by the
    // from-scratch computation in MarginalLatentParentReference.h (no memo, no incremental sums), with fresh COI and
    // strain-count distributions.
    namespace marginal = transmission_nets::model::transmission_process::marginal;
    auto latentTerm = [&](const std::shared_ptr<InfectionEvent>& infection, const std::vector<std::shared_ptr<InfectionEvent>>& observed) {
        const auto& memo = model->transmissionMemo;
        COIProbabilityImpl coiProb(std::make_shared<transmission_nets::core::parameters::Parameter<double>>(state->meanCOI->value()));
        const auto strainTable = TransmissionMemoImpl::Kernel::strainCountTable(state->meanStrainsTransmitted->value());
        ParentSetSizeLikelihoodImpl psp(std::make_shared<transmission_nets::core::parameters::Parameter<double>>(state->parentSetSizeProb->value()));
        const int n     = static_cast<int>(observed.size()) + 1;
        const auto variant = observed.empty() ? TransmissionMemoImpl::Variant::LatentParentOnly : TransmissionMemoImpl::Variant::ObservedAndLatentParents;
        marginal::LocusTable<MAX_COI, MAX_STRAINS> sums{};
        sums.fill(0.0);
        for (const auto& locus : infection->loci()) {
            const auto& child = infection->latentGenotype(locus)->value();
            std::vector<unsigned int> childAlleles;
            for (unsigned int a = 0; a < child.totalAlleles(); ++a) {
                if (child.allele(a)) {
                    childAlleles.push_back(a);
                }
            }
            if (childAlleles.empty()) {
                return -std::numeric_limits<double>::infinity();
            }
            const auto frequencies = state->alleleFrequencies->alleleFrequencies(locus)->value().frequencies();
            std::vector<const GeneticsImpl*> parents;
            for (const auto& p : observed) {
                parents.push_back(&p->latentGenotype(locus)->value());
            }
            const auto table = marginal::locusTable<MAX_COI, MAX_STRAINS, GeneticsImpl>(marginal::SourceOverlap(frequencies, childAlleles, MAX_COI), locus->totalAlleles(), [&](const GeneticsImpl& g) -> std::optional<std::array<double, MAX_STRAINS>> {
                const auto& term = memo->locusTerm(variant, child, std::span<const GeneticsImpl* const>(parents), &g);
                if (term.impossible) {
                    return std::nullopt;
                }
                return term.terms;
            });
            marginal::accumulate<MAX_COI, MAX_STRAINS>(sums, table);
        }
        return marginal::combine<MAX_COI, MAX_STRAINS>(sums, coiProb.value(), [&](const int m, const int np) { return TransmissionMemoImpl::Kernel::probNumStrains(strainTable, m, np); }, n, std::log(psp.value()(n)));
    };
    auto recompute = [&](const size_t j, ParentSetLliks* parentSets = nullptr) {
        const auto& infection = state->infections[j];
        const auto& latent    = state->latentParents[j];
        const auto& memo      = model->transmissionMemo;
        const auto& psp       = model->parentSetSizeLikelihood;
        const auto parents    = model->parentSets.at(infection->id())->value();
        const int totalNodes  = parents.size();

        std::vector<double> lliks{latentTerm(infection, {})};
        if (parentSets != nullptr) {
            parentSets->clear();
            parentSets->emplace_back(lliks.back(), transmission_nets::core::containers::ParentSet<InfectionEvent>{latent});
        }
        transmission_nets::core::utils::generators::CombinationIndicesGenerator comboGen;
        for (int i = 1; i <= MAX_PARENTS and i <= totalNodes; i++) {
            comboGen.reset(totalNodes, i);
            while (!comboGen.completed) {
                transmission_nets::core::containers::ParentSet<InfectionEvent> parentSet{};
                std::vector<std::shared_ptr<InfectionEvent>> observed;
                for (const auto& idx : comboGen.curr) {
                    parentSet.insert(parents.begin()[idx]);
                    observed.push_back(parents.begin()[idx]);
                }
                lliks.push_back(memo->observedParentsLogLikelihood(infection, parentSet) + std::log(psp->value()(parentSet.size())));
                if (parentSets != nullptr) {
                    parentSets->emplace_back(lliks.back(), parentSet);
                }
                lliks.push_back(latentTerm(infection, observed));
                if (parentSets != nullptr) {
                    parentSet.insert(latent);
                    parentSets->emplace_back(lliks.back(), parentSet);
                }
                comboGen.next();
            }
        }
        const double maxLlik = *std::max_element(lliks.begin(), lliks.end());
        if (maxLlik == -std::numeric_limits<double>::infinity()) {
            return maxLlik;
        }
        double sum = 0;
        for (const auto llik : lliks) {
            sum += std::exp(llik - maxLlik);
        }
        return maxLlik + std::log(sum);
    };

    // Marking each process dirty makes it recompute from whatever its cache holds, so stale entries surface
    // even when the process's own value was restored by a rejected proposal.
    int finite = 0;
    auto checkAll = [&](const std::string& when) {
        // The parent-set-size prior must follow its parameter; recompute() reads it from the model.
        ParentSetSizeLikelihoodImpl fresh(std::make_shared<transmission_nets::core::parameters::Parameter<double>>(state->parentSetSizeProb->value()));
        for (int k = 1; k <= MAX_PARENT_SET_SIZE; ++k) {
            EXPECT_DOUBLE_EQ(model->parentSetSizeLikelihood->value()(k), fresh.value()(k)) << when << ", parent set size " << k;
        }
        // The model's likelihood is the sum of its terms, each recomputed.
        const double modelLikelihood = model->getLikelihood();
        double termSum               = 0;
        for (const auto& observation : model->observationProcessLikelihoodList) {
            observation->setDirty();
            termSum += observation->value();
        }
        for (size_t j = 0; j < state->infections.size(); ++j) {
            // As the model reads it, which is stale if a change failed to mark the process dirty.
            const auto asRead = model->transmissionProcessList[j]->value();
            ParentSetLliks parentSets;
            const auto expected = recompute(j, &parentSets);
            // Parent sets are updated by the loci that changed rather than re-summed, so agreement is up to rounding.
            auto expectClose = [&](const double actual, const double reference, const std::string& how) {
                if (std::isfinite(reference)) {
                    EXPECT_NEAR(actual, reference, 1e-9 * std::max(1.0, std::abs(reference))) << when << ", node " << state->infections[j]->id() << ", " << how;
                } else {
                    EXPECT_EQ(actual, reference) << when << ", node " << state->infections[j]->id() << ", " << how;
                }
            };
            auto expectMatches = [&](const double actual, const std::string& how) { expectClose(actual, expected, how); };
            expectMatches(asRead, "as read");
            termSum += expected;
            model->transmissionProcessList[j]->setDirty();
            const auto cached = model->transmissionProcessList[j]->value();
            expectMatches(cached, "recomputed from its cache");
            // Every parent set's likelihood, as the parent-set logger reports it.
            const auto dist = model->transmissionProcessList[j]->calcParentSetDist();
            ASSERT_EQ(dist.parentSetLliks.size(), parentSets.size()) << when;
            for (size_t k = 0; k < parentSets.size(); ++k) {
                EXPECT_TRUE(dist.parentSetLliks[k].second == parentSets[k].second) << when << ", parent set " << k;
                expectClose(dist.parentSetLliks[k].first, parentSets[k].first, "parent set " + std::to_string(k));
            }
            if (std::isfinite(cached)) {
                ++finite;
            }
        }
        if (std::isfinite(termSum)) {
            EXPECT_NEAR(modelLikelihood, termSum, 1e-9 * std::max(1.0, std::abs(termSum))) << when << ", the model's likelihood";
        }
        // Marking the terms dirty above queued them in the model's accumulator; it must take them in before any restore
        // (an accumulator restored with targets still queued counts them twice).
        model->value();
        EXPECT_NEAR(model->getLikelihood(), modelLikelihood, 1e-9 * std::max(1.0, std::abs(modelLikelihood))) << when << ", the model's likelihood after re-evaluating its terms";
    };

    // Like a path proposal: one locus of a node and of one of its current parents change before a single evaluation.
    auto perturbWithParent = [&](const int node) {
        const auto& infection = state->infections[node];
        const auto parents    = model->parentSets.at(infection->id())->value();
        if (parents.empty()) {
            return;
        }
        const auto& parent = parents.begin()[node % parents.size()];
        const auto& locus  = infection->loci()[node % infection->loci().size()];
        auto childGenotype  = infection->latentGenotype(locus);
        auto parentGenotype = parent->latentGenotype(locus);
        auto step           = transmission_nets::core::samplers::Proposal(*model, *rng, perturbations);
        step.change(*parentGenotype, GeneticsImpl(randomBits(locus->totalAlleles())));
        step.change(*childGenotype, parentGenotype->value());
        if (std::isfinite(step.value()) && unif(*rng) < 0.5) {
            step.accept();
        } else {
            step.reject();
        }
        model->value();
    };

    // Like a whole-genome proposal: most loci of an infection change before a single evaluation, to those of a donor
    // so that the model mostly stays possible. Changing most loci makes parent sets be recomputed rather than updated.
    auto perturbGenome = [&](const auto& infection, const auto& donor) {
        const auto& loci = infection->loci();
        auto step        = transmission_nets::core::samplers::Proposal(*model, *rng, perturbations);
        // Sometimes the whole genome, so that the infection joins the donor's genome class.
        const double share = unif(*rng) < 0.5 ? 1.0 : 0.8;
        for (const auto& locus : loci) {
            if (unif(*rng) < share) {
                step.change(*infection->latentGenotype(locus), donor->latentGenotype(locus)->value());
            }
        }
        if (std::isfinite(step.value()) && unif(*rng) < 0.5) {
            step.accept();
        } else {
            step.reject();
        }
        model->value();
    };

    // Like a move of several infections at once: a few nodes each take another's genome, or most of it, before a
    // single evaluation, so that within one parent set classes split, merge and empty together.
    auto perturbGenomes = [&]() {
        boost::random::uniform_int_distribution<> pickCount(2, 4);
        std::vector<std::shared_ptr<InfectionEvent>> changed;
        auto step = transmission_nets::core::samplers::Proposal(*model, *rng, perturbations);
        for (int c = pickCount(*rng); c > 0; --c) {
            const auto& infection = state->infections[pickNode(*rng)];
            const auto& donor     = state->infections[pickNode(*rng)];
            if (std::ranges::find(changed, infection) != changed.end()) {
                continue;
            }
            changed.push_back(infection);
            const double share = unif(*rng) < 0.5 ? 1.0 : 0.8;
            for (const auto& locus : infection->loci()) {
                if (unif(*rng) < share) {
                    step.change(*infection->latentGenotype(locus), donor->latentGenotype(locus)->value());
                }
            }
        }
        if (std::isfinite(step.value()) && unif(*rng) < 0.5) {
            step.accept();
        } else {
            step.reject();
        }
        model->value();
    };

    for (int k = 0; k < 500; ++k) {
        if (k % 10 == 0) {
            checkAll("after " + std::to_string(k) + " perturbations");
        }
        if (k % 50 == 25) {
            // A real sweep, between checks, so that its moves meet the perturbations that follow.
            samplers.step();
        }
        const int node  = pickNode(*rng);
        const auto kind = k % 5;
        if (kind == 2) {
            auto duration = state->infections[node]->infectionDuration();
            perturb(duration, 1.0 + unif(*rng) * 58.0);
        } else if (kind == 3) {
            // Global parameters; allele frequencies and mean COI change only the parent sets with the latent parent.
            switch ((k / 5) % 4) {
                case 0:
                    perturb(state->parentSetSizeProb, 0.2 + 0.79 * unif(*rng));
                    break;
                case 1:
                    perturb(state->meanStrainsTransmitted, 1.0 + 3.0 * unif(*rng));
                    break;
                case 2:
                    perturb(state->meanCOI, 1.0 + 3.0 * unif(*rng));
                    break;
                default: {
                    auto locus = std::next(state->loci.begin(), static_cast<long>(k % state->loci.size()))->second;
                    std::vector<double> frequencies(locus->totalAlleles());
                    for (auto& frequency : frequencies) {
                        frequency = 0.05 + unif(*rng);
                    }
                    auto parameter = state->alleleFrequencies->alleleFrequencies(locus);
                    perturb(parameter, transmission_nets::core::datatypes::Simplex(frequencies));
                }
            }
            // Before another move dirties the processes for other reasons.
            checkAll("after global parameter change " + std::to_string(k));
        } else if (kind == 4) {
            perturbWithParent(node);
        } else if ((k / 5) % 3 == 1 and kind == 1) {
            perturbGenomes();
        } else if ((k / 5) % 3 == 2) {
            // A node takes most of another node's genome.
            perturbGenome(state->infections[node], state->infections[pickNode(*rng)]);
        } else {
            const auto& infection = state->infections[kind == 0 ? node : pickNode(*rng)];
            const auto& loci      = infection->loci();
            const auto& locus     = loci[k % loci.size()];
            auto genotype         = infection->latentGenotype(locus);
            perturb(genotype, GeneticsImpl(randomBits(locus->totalAlleles())));
        }
    }

    checkAll("at the end");

    // A rejected change to a node's own genotype is the one case nothing else repairs: its cache was cleared and
    // refilled under the rejected genotype, so the restore alone must bring the earlier entries back.
    for (size_t j = 0; j < state->infections.size(); ++j) {
        const auto& infection = state->infections[j];
        const auto& locus     = infection->loci()[j % infection->loci().size()];
        auto genotype         = infection->latentGenotype(locus);
        genotype->saveState(1);
        genotype->setValue(GeneticsImpl(randomBits(locus->totalAlleles())));
        model->value();
        genotype->restoreState(1);
        model->value();
    }
    checkAll("after rejected changes to each node's genotype");

    // Nested saved states: the node's own genotype is saved under one id and the allele frequencies at the same locus
    // under another, and both reach the node's process. In the second round the inner state is saved without a
    // change, so the process's inner cache is shared with the outer one while its own slot holds stale entries from
    // the first round. As in perturb(), an impossible state is restored at once.
    int nested = 0;
    for (int round = 0; round < 2; ++round) {
        for (size_t j = 0; j < state->infections.size(); ++j) {
            const auto& locus = state->infections[j]->loci()[(j + 1) % state->infections[j]->loci().size()];
            auto outer        = state->infections[j]->latentGenotype(locus);
            auto inner        = state->alleleFrequencies->alleleFrequencies(locus);
            const auto when   = "round " + std::to_string(round) + ", nested changes for node " + state->infections[j]->id();
            // The child takes another infection's genotype, which mostly keeps the model possible.
            outer->saveState(1);
            outer->setValue(state->infections[pickNode(*rng)]->latentGenotype(locus)->value());
            if (not std::isfinite(model->value())) {
                outer->restoreState(1);
                model->value();
                continue;
            }
            inner->saveState(2);
            const bool innerChanged = round == 0;
            if (innerChanged) {
                std::vector<double> frequencies(locus->totalAlleles());
                for (auto& frequency : frequencies) {
                    frequency = 0.05 + unif(*rng);
                }
                inner->setValue(transmission_nets::core::datatypes::Simplex(frequencies));
            }
            if (not std::isfinite(model->value())) {
                inner->restoreState(2);
                model->value();
            } else {
                ++nested;
                // Checking recomputes every process, which would give an unchanged inner state its own copy.
                if (innerChanged) {
                    checkAll("inside " + when);
                }
                if (j % 3 == 1) {
                    inner->restoreState(2);
                    model->value();
                }
            }
            if (j % 3 != 0) {
                outer->acceptState();
                inner->acceptState();
            } else {
                inner->restoreState(2);
                outer->restoreState(1);
            }
            model->value();
            checkAll("after " + when);
        }
    }
    EXPECT_GT(nested, static_cast<int>(state->infections.size()));
    EXPECT_GT(finite, 0);
}
}// namespace
