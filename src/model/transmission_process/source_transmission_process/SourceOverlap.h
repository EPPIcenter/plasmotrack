//
// The Source, at one locus: how a Latent parent's genotype G overlaps a child's genotype C.
//
// The Source draws a Latent parent as: COI k ~ P(k), then at each locus G = the distinct alleles of k i.i.d. draws from
// the locus' allele frequencies π. The transmission term depends on G only through U = G ∩ C and |G|, so what the
// model needs from the Source is P(G ∩ C = U, |G \ C| = w | k), which comes from the exponential generating function
// of k draws:
//   P(U, w | k) = k! [t^k x^w] prod_{i ∈ U} (e^{π_i t} - 1) prod_{i ∉ C} (1 + x (e^{π_i t} - 1)).
// See doc/marginal-latent-parent.md for the derivation.
//

#ifndef TRANSMISSION_NETWORKS_APP_SOURCEOVERLAP_H
#define TRANSMISSION_NETWORKS_APP_SOURCEOVERLAP_H

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include <limits>
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

}// namespace transmission_nets::model::transmission_process::marginal

#endif//TRANSMISSION_NETWORKS_APP_SOURCEOVERLAP_H
