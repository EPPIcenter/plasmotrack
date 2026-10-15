// Exploratory: would a chain prefer another chain's configuration of a few infections?
//
//   transplant INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP TARGET_CHECKPOINT SOURCE_CHECKPOINT ID [ID ...]
//
// Restores the hot chains of both checkpoints. For the listed infections, copies the source chain's durations (order),
// latent genotypes (genomes), or both into the target chain's state, and reports the change in the target chain's log
// posterior (and its prior and likelihood). Every change is undone after it is evaluated.

#include "core/io/parse_json.h"
#include "core/io/utils.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/EnsembleCheckpoint.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/State.h"
#include "impl/model/Model/config.h"

#include <fmt/core.h>

#include <string>
#include <vector>

using namespace transmission_nets;
using namespace transmission_nets::impl::Model;

int main(int argc, char** argv) {
    if (argc < 7) {
        fmt::print(stderr, "usage: {} INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP TARGET_CHECKPOINT SOURCE_CHECKPOINT ID [ID ...]\n", argv[0]);
        return 2;
    }
    const auto input = core::io::loadJSON(core::io::decompressGzipFile(argv[1]));
    const auto sIdp  = core::io::loadVectorFromFile<core::computation::Probability>(argv[2]);
    const auto aIdp  = core::io::loadVectorFromFile<core::computation::Probability>(argv[3]);
    auto rngA        = std::make_shared<checkpoint::Engine>(1);
    auto rngB        = std::make_shared<checkpoint::Engine>(2);
    auto target      = checkpoint::restore(checkpoint::load(argv[4]), rngA, input, sIdp, aIdp, false, 1);
    auto source      = checkpoint::restore(checkpoint::load(argv[5]), rngB, input, sIdp, aIdp, false, 1);
    const auto& t    = target->hotChain();
    const auto& s    = source->hotChain();

    auto index = [&](const State& state, const std::string& id) {
        for (size_t k = 0; k < state.infections.size(); ++k) {
            if (state.infections[k]->id() == id) {
                return k;
            }
        }
        throw std::out_of_range(id);
    };
    std::vector<size_t> chosen;
    for (int a = 6; a < argc; ++a) {
        chosen.push_back(index(*t.state, argv[a]));
    }

    const double base = t.model->value(), basePrior = t.model->getPrior(), baseLik = t.model->getLikelihood();
    fmt::print("target chain: log posterior {:.2f} (prior {:.2f}, likelihood {:.2f}); source chain: {:.2f}\n", base, basePrior, baseLik, s.model->value());
    for (const size_t k : chosen) {
        const auto& ti = t.state->infections[k];
        const auto& si = s.state->infections[index(*s.state, ti->id())];
        int differing  = 0;
        for (const auto& locus : ti->loci()) {
            const auto& a = ti->latentGenotype(locus)->value();
            const auto& b = si->latentGenotype(s.state->loci.at(locus->label))->value();
            for (unsigned int j = 0; j < a.totalAlleles(); ++j) {
                differing += a.allele(j) != b.allele(j);
            }
        }
        fmt::print("  {:>4}: infection time target {:7.1f} source {:7.1f}; genome alleles differing {}\n", ti->id(), ti->infectionTime(), si->infectionTime(), differing);
    }

    constexpr int stateId = 902;
    for (const auto& [label, order, genomes] : {std::tuple{"order only", true, false}, std::tuple{"genomes only", false, true}, std::tuple{"order and genomes", true, true}}) {
        for (const size_t k : chosen) {
            const auto& ti = t.state->infections[k];
            const auto& si = s.state->infections[index(*s.state, ti->id())];
            if (order) {
                ti->infectionDuration()->saveState(stateId);
                ti->infectionDuration()->setValue(si->infectionDuration()->value());
            }
            if (genomes) {
                for (const auto& locus : ti->loci()) {
                    ti->latentGenotype(locus)->saveState(stateId);
                    ti->latentGenotype(locus)->setValue(si->latentGenotype(s.state->loci.at(locus->label))->value());
                }
            }
        }
        const double v = t.model->value();
        fmt::print("  {:<18} change {:+9.2f} (prior {:+8.2f}, likelihood {:+9.2f})\n", label, v - base, t.model->getPrior() - basePrior, t.model->getLikelihood() - baseLik);
        for (const size_t k : chosen) {
            const auto& ti = t.state->infections[k];
            if (order) {
                ti->infectionDuration()->restoreState(stateId);
            }
            if (genomes) {
                for (const auto& locus : ti->loci()) {
                    ti->latentGenotype(locus)->restoreState(stateId);
                }
            }
        }
        t.model->value();
    }
    return 0;
}
