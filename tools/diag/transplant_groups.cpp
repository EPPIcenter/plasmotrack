// Exploratory: which part of another chain's state would a stuck chain prefer?
//
//   transplant_groups INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP TARGET_CHECKPOINT SOURCE_CHECKPOINT [TOP]
//
// Restores the hot chains of both checkpoints and copies groups of the source chain's state into the target chain's:
// every infection's genome, every duration, every false-positive and false-negative rate, every allele frequency, the
// global parameters, and everything together. Reports each change in the target's log posterior, then scans infections
// one at a time (genome; genome and duration; genome, duration and error rates) and lists the TOP largest gains.
// Every change is undone after it is evaluated.

#include "core/io/parse_json.h"
#include "core/io/utils.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/State.h"
#include "impl/model/Model/config.h"

#include <fmt/core.h>

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

using namespace transmission_nets;
using namespace transmission_nets::impl::Model;

namespace {
    constexpr int STATE_ID = 903;

    // A change that can be applied to the target and undone.
    struct Change {
        std::vector<std::function<void()>> apply, undo;
    };

    template<typename P, typename V>
    void add(Change& c, const std::shared_ptr<P>& param, const V& value) {
        c.apply.emplace_back([param, value] {
            param->saveState(STATE_ID);
            param->setValue(value);
        });
        c.undo.emplace_back([param] { param->restoreState(STATE_ID); });
    }
}// namespace

int main(int argc, char** argv) {
    if (argc < 6) {
        fmt::print(stderr, "usage: {} INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP TARGET_CHECKPOINT SOURCE_CHECKPOINT [TOP]\n", argv[0]);
        return 2;
    }
    const int top     = argc > 6 ? std::stoi(argv[6]) : 10;
    const auto input  = core::io::loadJSON(core::io::decompressGzipFile(argv[1]));
    const auto sIdp   = core::io::loadVectorFromFile<core::computation::Probability>(argv[2]);
    const auto aIdp   = core::io::loadVectorFromFile<core::computation::Probability>(argv[3]);
    auto target       = checkpoint::restore(checkpoint::load(argv[4]), std::make_shared<checkpoint::Engine>(1), input, sIdp, aIdp, false, 1);
    auto source       = checkpoint::restore(checkpoint::load(argv[5]), std::make_shared<checkpoint::Engine>(2), input, sIdp, aIdp, false, 1);
    const auto& t     = *target->hotChain().state;
    const auto& s     = *source->hotChain().state;
    auto& model       = *target->hotChain().model;
    const double base = model.value();
    fmt::print("target log posterior {:.2f} (prior {:.2f}, likelihood {:.2f}); source {:.2f}\n\n", base, model.getPrior(), model.getLikelihood(), source->hotChain().model->value());

    auto evaluate = [&](const Change& c) {
        const double prior0 = model.getPrior(), lik0 = model.getLikelihood();
        for (const auto& f : c.apply) {
            f();
        }
        const double v = model.value();
        const double dPrior = model.getPrior() - prior0, dLik = model.getLikelihood() - lik0;
        for (auto it = c.undo.rbegin(); it != c.undo.rend(); ++it) {
            (*it)();
        }
        model.value();
        return std::tuple{v - base, dPrior, dLik};
    };

    // Infections are in the same order in both states (both built from the same input).
    auto genome = [&](Change& c, size_t k) {
        const auto& ti = t.infections[k];
        const auto& si = s.infections[k];
        for (const auto& locus : ti->loci()) {
            const auto& value = si->latentGenotype(s.loci.at(locus->label))->value();
            if (value != ti->latentGenotype(locus)->value()) {
                add(c, ti->latentGenotype(locus), value);
            }
        }
    };
    auto duration = [&](Change& c, size_t k) { add(c, t.infections[k]->infectionDuration(), s.infections[k]->infectionDuration()->value()); };
    auto errors   = [&](Change& c, size_t k) {
        add(c, t.expectedFalsePositives[k], s.expectedFalsePositives[k]->value());
        add(c, t.expectedFalseNegatives[k], s.expectedFalseNegatives[k]->value());
    };
    auto frequencies = [&](Change& c) {
        for (const auto& [label, locus] : t.loci) {
            add(c, t.alleleFrequencies->alleleFrequencies(locus), s.alleleFrequencies->alleleFrequencies(s.loci.at(label))->value());
        }
    };
    auto globals = [&](Change& c) {
        add(c, t.meanCOI, s.meanCOI->value());
        add(c, t.meanStrainsTransmitted, s.meanStrainsTransmitted->value());
        add(c, t.parentSetSizeProb, s.parentSetSizeProb->value());
    };

    int differingGenomes = 0;
    for (size_t k = 0; k < t.infections.size(); ++k) {
        Change c;
        genome(c, k);
        differingGenomes += not c.apply.empty();
    }
    fmt::print("infections whose genomes differ between the chains: {} of {}\n", differingGenomes, t.infections.size());

    const std::vector<std::pair<std::string, std::function<void(Change&)>>> groups{
            {"all genomes", [&](Change& c) { for (size_t k = 0; k < t.infections.size(); ++k) genome(c, k); }},
            {"all durations", [&](Change& c) { for (size_t k = 0; k < t.infections.size(); ++k) duration(c, k); }},
            {"all error rates", [&](Change& c) { for (size_t k = 0; k < t.infections.size(); ++k) errors(c, k); }},
            {"allele frequencies", frequencies},
            {"global parameters", globals},
            {"genomes + durations", [&](Change& c) { for (size_t k = 0; k < t.infections.size(); ++k) { genome(c, k); duration(c, k); } }},
            {"everything", [&](Change& c) { for (size_t k = 0; k < t.infections.size(); ++k) { genome(c, k); duration(c, k); errors(c, k); } frequencies(c); globals(c); }},
    };
    for (const auto& [label, build] : groups) {
        Change c;
        build(c);
        const auto [d, dp, dl] = evaluate(c);
        fmt::print("  {:<22} change {:+9.2f} (prior {:+8.2f}, likelihood {:+9.2f})\n", label, d, dp, dl);
    }

    for (const auto& [label, parts] : std::vector<std::pair<std::string, int>>{{"genome", 1}, {"genome + duration", 2}, {"genome + duration + error rates", 3}}) {
        std::vector<std::tuple<double, double, double, size_t>> gains;
        for (size_t k = 0; k < t.infections.size(); ++k) {
            Change c;
            genome(c, k);
            if (parts >= 2) duration(c, k);
            if (parts >= 3) errors(c, k);
            const auto [d, dp, dl] = evaluate(c);
            gains.emplace_back(d, dp, dl, k);
        }
        std::ranges::sort(gains, std::greater{});
        fmt::print("\nlargest gains, one infection at a time ({}):\n", label);
        for (int i = 0; i < top and i < static_cast<int>(gains.size()); ++i) {
            const auto [d, dp, dl, k] = gains[i];
            fmt::print("  {:>8}  change {:+8.2f} (prior {:+7.2f}, likelihood {:+8.2f})  infection time target {:7.1f} source {:7.1f}\n", t.infections[k]->id(), d, dp, dl, t.infections[k]->infectionTime(), s.infections[k]->infectionTime());
        }
    }
    return 0;
}
