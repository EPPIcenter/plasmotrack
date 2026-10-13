//
// The expectation over a latent parent drawn from the Source of a parent set's transmission term, computed exactly.
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
//

#ifndef TRANSMISSION_NETWORKS_APP_MARGINALLATENTPARENT_H
#define TRANSMISSION_NETWORKS_APP_MARGINALLATENTPARENT_H

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace transmission_nets::model::transmission_process::marginal {

    // The Source's distribution, at one locus, of how a latent parent's genotype G overlaps a child's genotype C: for
    // each nonempty U ⊆ C (a bit mask over C's alleles in index order), each w = |G \ C| and each COI k,
    // P(G ∩ C = U, |G \ C| = w | k). U = ∅ is left out: such a latent parent makes the parent set impossible.
    class SourceOverlap {
    public:
        SourceOverlap() = default;
        SourceOverlap(const std::vector<double>& frequencies, const std::vector<unsigned int>& childAlleles, const int maxCoi)
            : childAlleles_(childAlleles), maxCoi_(maxCoi) {
            const auto numAlleles = static_cast<unsigned int>(frequencies.size());
            const auto childSize  = static_cast<unsigned int>(childAlleles.size());
            assert(childSize >= 1 and childSize < 32);
            const int K = maxCoi;
            maxW_       = std::min<unsigned int>(numAlleles - childSize, static_cast<unsigned int>(K));
            const size_t stride = K + 1;

            // 1 / b!, and k!.
            std::vector<double> invFactorial(stride, 1.0);
            factorial_.assign(stride, 1.0);
            for (int b = 1; b <= K; ++b) {
                factorial_[b]   = factorial_[b - 1] * b;
                invFactorial[b] = 1.0 / factorial_[b];
            }
            // Series of e^{πt} - 1 to degree K, into s.
            std::vector<double> s(stride);
            auto drawnSeries = [&](const double p) {
                s[0]         = 0.0;
                double power = 1.0;
                for (int b = 1; b <= K; ++b) {
                    power *= p;
                    s[b] = power * invFactorial[b];
                }
            };

            std::vector<bool> inChild(numAlleles, false);
            for (const auto a : childAlleles) {
                inChild[a] = true;
            }

            // R[w * stride + a]: alleles outside C, w of them drawn, coefficient of t^a.
            R_.assign((maxW_ + 1) * stride, 0.0);
            R_[0] = 1.0;
            for (unsigned int i = 0; i < numAlleles; ++i) {
                if (inChild[i] or frequencies[i] <= 0) {
                    continue;
                }
                drawnSeries(frequencies[i]);
                // R[w][a] is 0 for a < w (w distinct alleles need w draws), so only a >= w is updated, and only the
                // terms with a - b >= w - 1 are added; the zeros skipped would add nothing.
                for (int w = static_cast<int>(maxW_); w >= 1; --w) {
                    const double* below = R_.data() + (w - 1) * stride;
                    double* row         = R_.data() + w * stride;
                    for (int a = K; a >= w; --a) {
                        double add = 0;
                        for (int b = 1; b <= a - (w - 1); ++b) {
                            add += below[a - b] * s[b];
                        }
                        row[a] += add;
                    }
                }
            }

            childSeries_.resize(childSize * stride);
            for (unsigned int c = 0; c < childSize; ++c) {
                drawnSeries(frequencies[childAlleles[c]]);
                std::copy(s.begin(), s.end(), childSeries_.begin() + c * stride);
            }
            // A mask's probabilities are computed when first asked for: a parent set reads only the masks its observed
            // parents leave possible.
            offsets_.assign(1u << childSize, NOT_COMPUTED);
        }

        [[nodiscard]] double operator()(const unsigned int mask, const unsigned int w, const int k) const { return row(mask, w)[k]; }
        // The probabilities for k = 0..maxCoi, for one (U, w). Valid until the next call for a mask not asked for before.
        [[nodiscard]] const double* row(const unsigned int mask, const unsigned int w) const {
            if (offsets_[mask] == NOT_COMPUTED) {
                compute(mask);
            }
            return prob_.data() + offsets_[mask] + static_cast<size_t>(w) * (maxCoi_ + 1);
        }
        [[nodiscard]] unsigned int maxW() const { return maxW_; }
        [[nodiscard]] int maxCoi() const { return maxCoi_; }
        [[nodiscard]] const std::vector<unsigned int>& childAlleles() const { return childAlleles_; }

    private:
        static constexpr std::uint32_t NOT_COMPUTED = std::numeric_limits<std::uint32_t>::max();

        void compute(const unsigned int mask) const {
            const int K         = maxCoi_;
            const size_t stride = K + 1;
            const auto offset   = prob_.size();
            assert(offset < NOT_COMPUTED);
            offsets_[mask] = static_cast<std::uint32_t>(offset);
            prob_.resize(offset + (maxW_ + 1) * stride, 0.0);

            std::vector<double> E(stride, 0.0), next(stride);
            E[0] = 1.0;
            for (unsigned int c = 0; c < childAlleles_.size(); ++c) {
                if (not(mask >> c & 1u)) {
                    continue;
                }
                const double* s = childSeries_.data() + c * stride;
                std::fill(next.begin(), next.end(), 0.0);
                for (int a = 0; a <= K; ++a) {
                    if (E[a] == 0) {
                        continue;
                    }
                    for (int b = 1; a + b <= K; ++b) {
                        next[a + b] += E[a] * s[b];
                    }
                }
                std::swap(E, next);
            }
            // E[a] is 0 below |U| (each of U's alleles needs a draw) and row[k - a] below w, so a runs from |U|
            // to k - w; the zeros skipped would add nothing.
            const int u = std::popcount(mask);
            for (unsigned int w = 0; w <= maxW_; ++w) {
                const double* row = R_.data() + w * stride;
                double* out       = prob_.data() + offset + w * stride;
                for (int k = 1; k <= K; ++k) {
                    double coef = 0;
                    for (int a = u; a <= k - static_cast<int>(w); ++a) {
                        coef += E[a] * row[k - a];
                    }
                    out[k] = coef * factorial_[k];
                }
            }
        }

        std::vector<unsigned int> childAlleles_;
        int maxCoi_        = 0;
        unsigned int maxW_ = 0;
        std::vector<double> factorial_;
        std::vector<double> R_;
        std::vector<double> childSeries_;// e^{πt} - 1 to degree maxCoi for each child allele
        mutable std::vector<std::uint32_t> offsets_;// by mask: where its probabilities start in prob_
        mutable std::vector<double> prob_;
    };

    // log g_l(k, m) for one locus: entry k * MAX_STRAINS + (m - 1), for k = 0..MAX_COI (k = 0 is -inf).
    template<int MAX_COI, int MAX_STRAINS>
    using LocusTable = std::array<double, (MAX_COI + 1) * MAX_STRAINS>;

    // A representative latent-parent genotype with G ∩ C = U and |G \ C| = w: U's alleles and the first w alleles
    // outside C. The locus term depends on G only through U and |G|.
    template<typename Genetics>
    Genetics representative(const unsigned int totalAlleles, const std::vector<unsigned int>& childAlleles, const unsigned int mask, unsigned int w) {
        Genetics g(static_cast<int>(totalAlleles));
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
    template<int MAX_COI, int MAX_STRAINS, typename Genetics, typename LogTermFn>
    LocusTable<MAX_COI, MAX_STRAINS> locusTable(const SourceOverlap& overlap, const unsigned int totalAlleles, LogTermFn&& logTerm) {
        std::array<double, (MAX_COI + 1) * MAX_STRAINS> g{};
        g.fill(0.0);
        const auto& childAlleles = overlap.childAlleles();
        const unsigned int numMasks = 1u << childAlleles.size();
        for (unsigned int mask = 1; mask < numMasks; ++mask) {
            for (unsigned int w = 0; w <= overlap.maxW(); ++w) {
                const auto r = static_cast<int>(std::popcount(mask) + w);
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
                const std::optional<std::array<double, MAX_STRAINS>> term = logTerm(representative<Genetics>(totalAlleles, childAlleles, mask, w));
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

    // Adds a locus' table to a running sum (-inf absorbs).
    template<int MAX_COI, int MAX_STRAINS>
    void accumulate(LocusTable<MAX_COI, MAX_STRAINS>& sums, const LocusTable<MAX_COI, MAX_STRAINS>& table) {
        for (size_t i = 0; i < sums.size(); ++i) {
            sums[i] += table[i];
        }
    }

}// namespace transmission_nets::model::transmission_process::marginal

#endif//TRANSMISSION_NETWORKS_APP_MARGINALLATENTPARENT_H
