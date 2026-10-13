// Exploratory: how well would HMC sample one locus' allele frequencies, compared with the samplers we have?
//
//   frequency_hmc INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP CHECKPOINT LOCUS METHOD BUDGET SEED OUT_CSV
//   frequency_hmc INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP CHECKPOINT --list
//
// Restores the checkpoint's hot chain and samples the conditional distribution of LOCUS' frequencies with everything
// else held fixed. METHOD is
//   salt    SALTSampler as the model registers it (k proposals per update)
//   walk    SimplexLogRatioSampler (one proposal per update)
//   hmcL    HMC with L leapfrog steps (e.g. hmc10), in additive log-ratio coordinates, with a diagonal mass matrix and
//           step size (dual averaging to acceptance 0.8) tuned during warm-up; gradients by central differences of the
//           model's log posterior.
// BUDGET counts target evaluations for salt and walk, and gradients for HMC. The first 25% of it is warm-up
// (adaptation); the rest samples. OUT_CSV has a row per update: phase, evaluations, gradients, then the frequencies.

#include "core/io/parse_json.h"
#include "core/io/utils.h"
#include "core/samplers/general/SALTSampler.h"
#include "core/samplers/general/SimplexLogRatioSampler.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/State.h"
#include "impl/model/Model/config.h"

#include <boost/random.hpp>
#include <fmt/core.h>
#include <fmt/os.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace transmission_nets;
using namespace transmission_nets::impl::Model;

namespace {
    constexpr int STATE_ID = 904;

    // The model's log posterior, counting evaluations that recompute it (not those returning the cached value).
    struct CountingTarget {
        Model* model;
        long evaluations = 0;
        double value() {
            evaluations += model->isDirty();
            return model->value();
        }
        [[nodiscard]] bool isDirty() const { return model->isDirty(); }
    };

    using Simplex = core::datatypes::Simplex;
    using SimplexParameter = core::parameters::Parameter<Simplex>;

    struct Hmc {
        std::shared_ptr<SimplexParameter> freqs;
        CountingTarget* target;
        std::shared_ptr<checkpoint::Engine> rng;
        int leapfrogSteps;
        size_t k;
        long gradients = 0;
        boost::random::normal_distribution<> normal{0, 1};
        boost::random::uniform_01<> uniform{};

        // z_i = log(π_i / π_k), i < k.
        std::vector<double> toZ(const std::vector<double>& pi) const {
            std::vector<double> z(k - 1);
            for (size_t i = 0; i + 1 < k; ++i) {
                z[i] = std::log(pi[i]) - std::log(pi[k - 1]);
            }
            return z;
        }
        std::vector<double> toPi(const std::vector<double>& z) const {
            double maxZ = 0;
            for (const double v : z) {
                maxZ = std::max(maxZ, v);
            }
            std::vector<double> pi(k);
            double sum = std::exp(-maxZ);
            pi[k - 1]  = sum;
            for (size_t i = 0; i + 1 < k; ++i) {
                pi[i] = std::exp(z[i] - maxZ);
                sum += pi[i];
            }
            for (auto& p : pi) {
                p /= sum;
            }
            return pi;
        }
        // Log density of z: the model's log posterior at π(z) plus log |dπ/dz| = Σ log π_i.
        double logDensity(const std::vector<double>& z) {
            const auto pi = toPi(z);
            double jac    = 0;
            for (const double p : pi) {
                if (p < MIN_ALLELE_FREQUENCY) {
                    return -std::numeric_limits<double>::infinity();
                }
                jac += std::log(p);
            }
            freqs->saveState(STATE_ID);
            freqs->setValue(Simplex(pi));
            const double v = target->value();
            freqs->acceptState();
            return v + jac;
        }
        bool gradient(const std::vector<double>& z, std::vector<double>& g) {
            ++gradients;
            constexpr double h = 1e-5;
            auto w             = z;
            for (size_t i = 0; i + 1 < k; ++i) {
                w[i]            = z[i] + h;
                const double up = logDensity(w);
                w[i]            = z[i] - h;
                const double dn = logDensity(w);
                w[i]            = z[i];
                if (not std::isfinite(up) or not std::isfinite(dn)) {
                    return false;
                }
                g[i] = (up - dn) / (2 * h);
            }
            return true;
        }
    };
}// namespace

int main(int argc, char** argv) {
    if (argc < 6) {
        fmt::print(stderr, "usage: {} INPUT SYMPTOMATIC_IDP ASYMPTOMATIC_IDP CHECKPOINT (--list | LOCUS METHOD BUDGET SEED OUT_CSV)\n", argv[0]);
        return 2;
    }
    const auto input  = core::io::loadJSON(core::io::decompressGzipFile(argv[1]));
    const auto sIdp   = core::io::loadVectorFromFile<core::computation::Probability>(argv[2]);
    const auto aIdp   = core::io::loadVectorFromFile<core::computation::Probability>(argv[3]);
    auto ensemble     = checkpoint::restore(checkpoint::load(argv[4]), std::make_shared<checkpoint::Engine>(1), input, sIdp, aIdp, false, 1);
    auto& chain       = ensemble->hotChain();
    auto& state       = *chain.state;
    auto& model       = *chain.model;
    if (std::string(argv[5]) == "--list") {
        for (const auto& [label, locus] : state.loci) {
            const auto f = state.alleleFrequencies->alleleFrequencies(locus)->value().frequencies();
            fmt::print("{}\t{}\t{:.4f}\n", label, f.size(), *std::ranges::min_element(f));
        }
        return 0;
    }
    if (argc < 10) {
        fmt::print(stderr, "missing arguments\n");
        return 2;
    }
    const std::string label = argv[5], method = argv[6];
    const long budget       = std::stol(argv[7]);
    auto rng                = std::make_shared<checkpoint::Engine>(std::stoul(argv[8]));
    auto out                = fmt::output_file(argv[9]);
    auto freqs              = state.alleleFrequencies->alleleFrequencies(state.loci.at(label));
    const size_t k          = freqs->value().totalElements();
    auto target             = std::make_shared<CountingTarget>(CountingTarget{&model});
    const long warmup       = budget / 4;
    model.value();

    out.print("phase,evaluations,gradients");
    for (size_t i = 0; i < k; ++i) {
        out.print(",p{}", i);
    }
    out.print("\n");
    auto record = [&](const bool sampling, const long gradients) {
        out.print("{},{},{}", sampling ? "sampling" : "warmup", target->evaluations, gradients);
        for (const double p : freqs->value().frequencies()) {
            out.print(",{:.8g}", p);
        }
        out.print("\n");
    };

    const auto start = std::chrono::steady_clock::now();
    if (method == "salt" or method == "walk") {
        std::unique_ptr<core::samplers::AbstractSampler> sampler;
        if (method == "salt") {
            auto s = std::make_unique<core::samplers::SALTSampler<CountingTarget, checkpoint::Engine>>(freqs, target, rng, 1, .1, 2);
            s->setLowerLimit_(MIN_ALLELE_FREQUENCY);
            sampler = std::move(s);
        } else {
            sampler = std::make_unique<core::samplers::general::SimplexLogRatioSampler<CountingTarget, checkpoint::Engine>>(freqs, target, rng, MIN_ALLELE_FREQUENCY);
        }
        while (target->evaluations < budget) {
            const bool sampling = target->evaluations >= warmup;
            sampler->update();
            if (not sampling) {
                sampler->adapt();
            }
            record(sampling, 0);
        }
    } else if (method.starts_with("hmc")) {
        Hmc hmc{.freqs = freqs, .target = target.get(), .rng = rng, .leapfrogSteps = std::stoi(method.substr(3)), .k = k};
        const size_t d = k - 1;
        auto z         = hmc.toZ(freqs->value().frequencies());
        std::vector<double> g(d), invMass(d, 1.0);
        double logP = hmc.logDensity(z);
        if (not hmc.gradient(z, g)) {
            fmt::print(stderr, "no gradient at the start\n");
            return 1;
        }
        // Dual averaging (Hoffman & Gelman 2014) towards acceptance 0.8.
        double logEps = std::log(0.1), logEpsBar = 0, hBar = 0;
        const double mu = std::log(10 * 0.1);
        long iteration  = 0;
        // Welford variance of z over the second half of warm-up, for the mass matrix.
        std::vector<double> zMean(d, 0.0), zM2(d, 0.0);
        long zCount        = 0;
        bool massSet       = false;
        long accepted      = 0, proposed = 0;
        while (hmc.gradients < budget) {
            const bool sampling = hmc.gradients >= warmup;
            ++iteration;
            double eps = std::exp(sampling ? logEpsBar : logEps) * (0.8 + 0.4 * hmc.uniform(*rng));
            std::vector<double> p(d), zNew = z, gNew = g;
            double kinetic = 0;
            for (size_t i = 0; i < d; ++i) {
                p[i] = hmc.normal(*rng) / std::sqrt(invMass[i]);
                kinetic += 0.5 * p[i] * p[i] * invMass[i];
            }
            const double h0 = -logP + kinetic;
            bool ok         = true;
            double logPNew  = logP;
            for (int step = 0; step < hmc.leapfrogSteps and ok; ++step) {
                for (size_t i = 0; i < d; ++i) {
                    p[i] += 0.5 * eps * gNew[i];
                    zNew[i] += eps * invMass[i] * p[i];
                }
                logPNew = hmc.logDensity(zNew);
                ok      = std::isfinite(logPNew) and hmc.gradient(zNew, gNew);
                if (not std::isfinite(logPNew)) {
                    ++hmc.gradients;// a step that left the support still costs the budget a step
                }
                if (ok) {
                    for (size_t i = 0; i < d; ++i) {
                        p[i] += 0.5 * eps * gNew[i];
                    }
                }
            }
            double acceptProb = 0;
            if (ok) {
                double kineticNew = 0;
                for (size_t i = 0; i < d; ++i) {
                    kineticNew += 0.5 * p[i] * p[i] * invMass[i];
                }
                acceptProb = std::min(1.0, std::exp(h0 - (-logPNew + kineticNew)));
            }
            ++proposed;
            if (hmc.uniform(*rng) < acceptProb) {
                ++accepted;
                z    = zNew;
                g    = gNew;
                logP = logPNew;
            }
            // Leave the parameter at the chain's current point.
            logP = hmc.logDensity(z);
            if (not sampling) {
                const double eta = 1.0 / (iteration + 10);
                hBar             = (1 - eta) * hBar + eta * (0.8 - acceptProb);
                logEps           = mu - std::sqrt(static_cast<double>(iteration)) / 0.05 * hBar;
                const double w   = std::pow(static_cast<double>(iteration), -0.75);
                logEpsBar        = w * logEps + (1 - w) * logEpsBar;
                if (hmc.gradients >= warmup / 2) {
                    ++zCount;
                    for (size_t i = 0; i < d; ++i) {
                        const double delta = z[i] - zMean[i];
                        zMean[i] += delta / zCount;
                        zM2[i] += delta * (z[i] - zMean[i]);
                    }
                }
                if (not massSet and hmc.gradients >= warmup * 3 / 4 and zCount > 10) {
                    for (size_t i = 0; i < d; ++i) {
                        invMass[i] = std::max(zM2[i] / (zCount - 1), 1e-8);
                    }
                    massSet = true;
                    // Restart the step-size adaptation for the new metric.
                    iteration = 0;
                    hBar      = 0;
                    logEpsBar = 0;
                }
            }
            record(sampling, hmc.gradients);
        }
        fmt::print(stderr, "hmc acceptance {:.3f}, step {:.4g}, gradients {}, evaluations {}\n",
                   static_cast<double>(accepted) / proposed, std::exp(logEpsBar), hmc.gradients, target->evaluations);
    } else {
        fmt::print(stderr, "unknown method {}\n", method);
        return 2;
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    fmt::print(stderr, "{} on {} (k = {}): {} evaluations in {:.1f} s ({:.2f} ms each)\n", method, label, k, target->evaluations, seconds, 1000 * seconds / target->evaluations);
    return 0;
}
