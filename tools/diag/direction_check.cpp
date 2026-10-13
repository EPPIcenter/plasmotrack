// Exploratory: which direction of the edge between two infections does the model favour?
//
//   direction_check INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP CHECKPOINT A B
//
// Restores the hot chain of CHECKPOINT. With latent parents summed out of the model, an infection x's transmission
// term is log( W0 + sum_S E_{L~Source}[ W(S ∪ {L}) ] ), which depends on the order but on no latent-parent genome. The
// tool reports the pair's terms and the duration prior in the chain's order, then with one of the two moved to the
// other side of the other, at several durations. Every change is undone after it is evaluated.

#include "core/io/parse_json.h"
#include "core/io/utils.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/State.h"
#include "impl/model/Model/config.h"

#include <boost/random.hpp>
#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace transmission_nets;
using namespace transmission_nets::impl::Model;

namespace {
    constexpr int STATE_ID = 901;


    double logSumExp(const std::vector<double>& xs) {
        double m = -std::numeric_limits<double>::infinity();
        for (const double x : xs) {
            m = std::max(m, x);
        }
        if (m == -std::numeric_limits<double>::infinity()) {
            return m;
        }
        double s = 0;
        for (const double x : xs) {
            s += std::exp(x - m);
        }
        return m + std::log(s);
    }

    struct Checker {
        std::shared_ptr<State> state;
        std::shared_ptr<Model> model;

        // Each parent set's term, keyed by its sampled parents ("" for the latent parent alone, "+L" appended to
        // those with it).
        std::map<std::string, double> sets_{};

        // x's transmission term in the current order. Fills sets_.
        std::pair<double, double> marginal(const size_t x) {
            model->value();
            const auto dist   = model->transmissionProcessList[x]->calcParentSetDist();
            const auto latent = state->latentParents[x];
            sets_.clear();
            for (const auto& [llik, ps] : dist.parentSetLliks) {
                std::vector<std::string> members;
                for (const auto& m : ps) {
                    if (m != latent) {
                        members.push_back(m->id());
                    }
                }
                std::ranges::sort(members);
                std::string key;
                for (const auto& m : members) {
                    key += (key.empty() ? "" : "+") + m;
                }
                if (ps.contains(latent) and not key.empty()) {
                    key += "+L";
                }
                sets_[key] = llik;
            }
            return {model->transmissionProcessList[x]->value(), std::numeric_limits<double>::infinity()};
        }

        // From sets_: P(no sampled parent), and the most probable sampled parents with their probabilities.
        std::string summary() const {
            std::vector<double> all;
            for (const auto& [key, v] : sets_) {
                all.push_back(v);
            }
            const double total = logSumExp(all);
            std::map<std::string, double> parentProb;
            double none = 0;
            for (const auto& [key, v] : sets_) {
                const double p = std::exp(v - total);
                if (key.empty()) {
                    none += p;
                    continue;
                }
                std::string rest = key.ends_with("+L") ? key.substr(0, key.size() - 2) : key;
                size_t start = 0;
                while (start <= rest.size()) {
                    const size_t plus = rest.find('+', start);
                    parentProb[rest.substr(start, plus == std::string::npos ? std::string::npos : plus - start)] += p;
                    if (plus == std::string::npos) {
                        break;
                    }
                    start = plus + 1;
                }
            }
            std::vector<std::pair<double, std::string>> ranked;
            for (const auto& [k, p] : parentProb) {
                ranked.emplace_back(p, k);
            }
            std::ranges::sort(ranked, std::greater{});
            std::string out = fmt::format("P(no sampled parent) {:.2f}", none);
            for (size_t i = 0; i < ranked.size() and i < 3; ++i) {
                out += fmt::format(", {} {:.2f}", ranked[i].second, ranked[i].first);
            }
            return out;
        }
    };
}// namespace

int main(int argc, char** argv) {
    if (argc < 7) {
        fmt::print(stderr, "usage: {} INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP CHECKPOINT A B\n", argv[0]);
        return 2;
    }
    const auto input = core::io::loadJSON(core::io::decompressGzipFile(argv[1]));
    const auto sIdp  = core::io::loadVectorFromFile<core::computation::Probability>(argv[2]);
    const auto aIdp  = core::io::loadVectorFromFile<core::computation::Probability>(argv[3]);
    const auto saved = checkpoint::load(argv[4]);
    auto runRng      = std::make_shared<checkpoint::Engine>(1);
    auto ensemble    = checkpoint::restore(saved, runRng, input, sIdp, aIdp, false, 1);
    const auto& hot  = ensemble->hotChain();
    Checker c{hot.state, hot.model};

    auto index = [&](const std::string& id) {
        for (size_t k = 0; k < c.state->infections.size(); ++k) {
            if (c.state->infections[k]->id() == id) {
                return k;
            }
        }
        throw std::out_of_range(id);
    };
    const size_t a = index(argv[5]), b = index(argv[6]);
    const auto& A = c.state->infections[a];
    const auto& B = c.state->infections[b];

    auto evaluate = [&](const char* label) {
        const auto [ma, essA] = c.marginal(a);
        const auto sumA       = c.summary();
        const auto [mb, essB] = c.marginal(b);
        const auto sumB       = c.summary();
        const double prior    = c.model->getPrior();
        fmt::print("  {:<44} {} first. marginal {} {:9.2f}   {} {:9.2f}   sum {:9.2f}  prior {:9.2f}  sum + prior {:9.2f}\n", label,
                   A->infectionTime() < B->infectionTime() ? A->id() : B->id(), A->id(), ma, B->id(), mb, ma + mb, prior, ma + mb + prior);
        fmt::print("      averaged over L:  {}: {}\n                        {}: {}\n", A->id(), sumA, B->id(), sumB);
        return ma + mb + prior;
    };

    fmt::print("A = {}, B = {}; infection times A {:.1f}, B {:.1f}\n", A->id(), B->id(), A->infectionTime(), B->infectionTime());
    const double now = evaluate("as in the chain");
    for (const auto& [i, j] : {std::pair{a, b}, std::pair{b, a}}) {
        const auto& I     = c.state->infections[i];
        const auto& J     = c.state->infections[j];
        const bool iFirst = I->infectionTime() < J->infectionTime();
        const double maxDuration = static_cast<double>(I->isSymptomatic() ? sIdp.size() : aIdp.size());
        const double boundary    = I->observationTime()->value() - J->infectionTime();
        const double lo = iFirst ? 1.0 : std::max(1.0, boundary), hi = iFirst ? std::min(maxDuration, boundary) : maxDuration;
        if (hi <= lo) {
            continue;
        }
        for (const double q : {0.25, 0.5, 0.75}) {
            const double duration = lo + q * (hi - lo);
            I->infectionDuration()->saveState(STATE_ID);
            I->infectionDuration()->setValue(duration);
            c.model->value();
            if (std::isfinite(c.model->getPrior())) {
                const double v = evaluate(fmt::format("{} moved {} {} (duration {:.1f})", I->id(), iFirst ? "after" : "before", J->id(), duration).c_str());
                fmt::print("    change from the chain's order: {:+.2f}\n", v - now);
            } else {
                fmt::print("  {} moved {} {} (duration {:.1f}): duration prior is zero\n", I->id(), iFirst ? "after" : "before", J->id(), duration);
            }
            I->infectionDuration()->restoreState(STATE_ID);
            c.model->value();
        }
    }
    return 0;
}
