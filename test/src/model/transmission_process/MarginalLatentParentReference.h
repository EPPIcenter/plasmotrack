//
// The expectation over a latent parent drawn from the Source of a parent set's transmission term, computed exactly
// and plainly: the reference the tests check the Transmission memo's latent locus tables against.
// See doc/marginal-latent-parent.md for the derivation.
//
// For a child with latent genotypes C_l and a parent set S plus the latent parent L, the transmission term is
//   P(C | S ∪ {L}) = sum_m P(m | n) prod_l f_l(m, G_l)        (n = |S| + 1, m = strains transmitted)
// and the Source draws L as: COI k ~ P(k), then at each locus G_l = the distinct alleles of k draws from the locus'
// allele frequencies. f_l depends on G_l only through U = G_l ∩ C_l and r = |G_l|, and loci are independent given k,
// so
//   E_L[ P(C | S ∪ {L}) ] = sum_k P(k) sum_m P(m | n) prod_l g_l(k, m),
//   g_l(k, m) = sum_{U ⊆ C_l, U ≠ ∅} sum_w P(G ∩ C_l = U, |G \ C_l| = w | k) f_l(m; U, |U| + w).
// P(U, w | k) comes from the exponential generating function of k draws:
//   P(U, w | k) = k! [t^k x^w] prod_{i ∈ U} (e^{π_i t} - 1) prod_{i ∉ C} (1 + x (e^{π_i t} - 1)).
// With a transmission leak ε > 0, U = ∅ is possible for the latent parent alone, so the sum over U includes it.
//
// leakyLocusLogTerm() is the locus term with a leak, written plainly: the pool over every allele, then the strains
// drawn one at a time, tracking the set of alleles drawn so far.
//

#ifndef TRANSMISSION_NETWORKS_TEST_MARGINALLATENTPARENTREFERENCE_H
#define TRANSMISSION_NETWORKS_TEST_MARGINALLATENTPARENTREFERENCE_H

#include "core/datatypes/Genotype.h"
#include "model/transmission_process/source_transmission_process/SourceOverlap.h"

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace transmission_nets::model::transmission_process::marginal {

    // log g_l(k, m) for one locus: entry k * MAX_STRAINS + (m - 1), for k = 0..MAX_COI (k = 0 is -inf).
    template<int MAX_COI, int MAX_STRAINS>
    using LocusTable = std::array<double, (MAX_COI + 1) * MAX_STRAINS>;

    // A representative latent-parent genotype with G ∩ C = U and |G \ C| = w: U's alleles and the first w alleles
    // outside C. The locus term depends on G only through U and |G|.
    inline core::datatypes::Genotype representative(const unsigned int totalAlleles, const std::vector<unsigned int>& childAlleles, const unsigned int mask, unsigned int w) {
        core::datatypes::Genotype g(static_cast<int>(totalAlleles));
        g.reset();// the size constructor sets the last allele
        std::bitset<256> inChild;
        for (unsigned int c = 0; c < childAlleles.size(); ++c) {
            inChild.set(childAlleles[c]);
            if (mask >> c & 1u) {
                g.set(childAlleles[c]);
            }
        }
        for (unsigned int i = 0; i < totalAlleles and w > 0; ++i) {
            if (not inChild[i]) {
                g.set(i);
                --w;
            }
        }
        return g;
    }

    // log g(k, m) at one locus. `logTerm(latentGenotype)` gives the model's locus term for the parent set with the
    // latent parent carrying that genotype: log f(m) for m = 1..MAX_STRAINS, or nullopt when the parent set is
    // impossible.
    template<int MAX_COI, int MAX_STRAINS, typename LogTermFn>
    LocusTable<MAX_COI, MAX_STRAINS> locusTable(const SourceOverlap& overlap, const unsigned int totalAlleles, LogTermFn&& logTerm) {
        std::array<double, (MAX_COI + 1) * MAX_STRAINS> g{};
        g.fill(0.0);
        const auto& childAlleles = overlap.childAlleles();
        const unsigned int numMasks = 1u << childAlleles.size();
        // U = ∅ only contributes with a leak and the latent parent alone; logTerm says it is impossible otherwise.
        for (unsigned int mask = 0; mask < numMasks; ++mask) {
            for (unsigned int w = 0; w <= overlap.maxW(); ++w) {
                const auto r = static_cast<int>(std::popcount(mask) + w);
                if (r == 0) {
                    continue;// a latent parent carries at least one allele
                }
                if (r > MAX_COI) {
                    break;
                }
                // Any probability mass at all for this (U, w)?
                bool any = false;
                for (int k = r; k <= MAX_COI and not any; ++k) {
                    any = overlap(mask, w, k) > 0;
                }
                if (not any) {
                    continue;
                }
                const std::optional<std::array<double, MAX_STRAINS>> term = logTerm(representative(totalAlleles, childAlleles, mask, w));
                if (not term) {
                    continue;
                }
                std::array<double, MAX_STRAINS> f{};
                for (int m = 0; m < MAX_STRAINS; ++m) {
                    f[m] = std::exp((*term)[m]);
                }
                for (int k = r; k <= MAX_COI; ++k) {
                    const double p = overlap(mask, w, k);
                    if (p == 0) {
                        continue;
                    }
                    for (int m = 0; m < MAX_STRAINS; ++m) {
                        g[k * MAX_STRAINS + m] += p * f[m];
                    }
                }
            }
        }
        LocusTable<MAX_COI, MAX_STRAINS> out{};
        for (size_t i = 0; i < g.size(); ++i) {
            out[i] = g[i] > 0 ? std::log(g[i]) : -std::numeric_limits<double>::infinity();
        }
        return out;
    }

    // log E_L[ P(C | S ∪ {L}) ] including the parent-set-size prior, from the per-locus tables summed over the child's
    // loci (logSums, -inf where any locus is zero). coiLogProbs[k] = log P(k); logStrains(m, n) = log P(m | n).
    template<int MAX_COI, int MAX_STRAINS, typename StrainFn>
    double combine(const LocusTable<MAX_COI, MAX_STRAINS>& logSums, const std::array<double, MAX_COI + 1>& coiLogProbs, StrainFn&& logStrains, const int numParents, const double logParentSetSizeProb) {
        double maxTerm = -std::numeric_limits<double>::infinity();
        std::array<double, (MAX_COI + 1) * MAX_STRAINS> terms{};
        terms.fill(-std::numeric_limits<double>::infinity());
        for (int k = 1; k <= MAX_COI; ++k) {
            for (int m = numParents; m <= MAX_STRAINS; ++m) {
                const double s = logSums[k * MAX_STRAINS + m - 1];
                if (s == -std::numeric_limits<double>::infinity() or coiLogProbs[k] == -std::numeric_limits<double>::infinity()) {
                    continue;
                }
                const double t                  = coiLogProbs[k] + logStrains(m, numParents) + s;
                terms[k * MAX_STRAINS + m - 1]  = t;
                maxTerm                         = std::max(maxTerm, t);
            }
        }
        if (maxTerm == -std::numeric_limits<double>::infinity()) {
            return maxTerm;
        }
        double sum = 0;
        for (const double t : terms) {
            sum += std::exp(t - maxTerm);
        }
        return maxTerm + std::log(sum) + logParentSetSizeProb;
    }

    // log f(m), m = 1..MAX_STRAINS, at one locus with a transmission leak: m strains drawn from the pool
    //   p_j = (1 - leak) (1/n) sum_{q in parents} 1[j in G_q] / |G_q| + leak π_j
    // show exactly the child's alleles. `parents` are every parent of the set, the latent parent included (last when
    // `latentWithObserved`, which applies the rule that it must share an allele with the child). nullopt when no m is
    // possible. Computed by following the set of child alleles drawn so far, strain by strain; a draw outside the child's
    // alleles ends the path.
    template<int MAX_STRAINS>
    std::optional<std::array<double, MAX_STRAINS>> leakyLocusLogTerm(const core::datatypes::Genotype& child, const std::vector<const core::datatypes::Genotype*>& parents, const std::vector<double>& frequencies,
                                                                     const double leak, const bool latentWithObserved) {
        const unsigned int A = child.totalAlleles();
        std::vector<unsigned int> childAlleles;
        for (unsigned int a = 0; a < A; ++a) {
            if (child.allele(a)) {
                childAlleles.push_back(a);
            }
        }
        if (childAlleles.empty()) {
            return std::nullopt;
        }
        if (latentWithObserved and core::datatypes::Genotype::truePositiveCount(*parents.back(), child) == 0) {
            return std::nullopt;
        }
        std::vector<double> pool(A, 0.0);
        for (const auto* parent : parents) {
            const double size = parent->totalPositiveCount();
            for (unsigned int a = 0; a < A; ++a) {
                if (parent->allele(a)) {
                    pool[a] += 1.0 / size / static_cast<double>(parents.size());
                }
            }
        }
        for (unsigned int a = 0; a < A; ++a) {
            pool[a] = (1 - leak) * pool[a] + leak * frequencies[a];
        }
        const unsigned int c    = childAlleles.size();
        const unsigned int full = (1u << c) - 1;
        std::vector<double> drawn(1u << c, 0.0), next(1u << c);
        drawn[0] = 1.0;
        std::array<double, MAX_STRAINS> out{};
        bool any = false;
        for (int m = 1; m <= MAX_STRAINS; ++m) {
            std::fill(next.begin(), next.end(), 0.0);
            for (unsigned int set = 0; set <= full; ++set) {
                for (unsigned int i = 0; i < c; ++i) {
                    next[set | 1u << i] += drawn[set] * pool[childAlleles[i]];
                }
            }
            std::swap(drawn, next);
            out[m - 1] = drawn[full] > 0 ? std::log(drawn[full]) : -std::numeric_limits<double>::infinity();
            any        = any or drawn[full] > 0;
        }
        if (not any) {
            return std::nullopt;
        }
        return out;
    }

    // Adds a locus' table to a running sum (-inf absorbs).
    template<int MAX_COI, int MAX_STRAINS>
    void accumulate(LocusTable<MAX_COI, MAX_STRAINS>& sums, const LocusTable<MAX_COI, MAX_STRAINS>& table) {
        for (size_t i = 0; i < sums.size(); ++i) {
            sums[i] += table[i];
        }
    }

}// namespace transmission_nets::model::transmission_process::marginal

#endif//TRANSMISSION_NETWORKS_TEST_MARGINALLATENTPARENTREFERENCE_H
