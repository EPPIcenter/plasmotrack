//
// The Transmission kernel: the arithmetic of the Transmission process for one Locus or one Parent set, a pure function
// of genotypes, Allele frequencies and the count distributions. It holds no Infection, no Parameter and no cache; the
// Transmission memo keeps its results, and each Infection's Transmission term sums them.
//
// For a Parent set S of a child with genotype C_l at locus l, the term is sum_m P(m | n) prod_l f_l(m), n = |S| and m
// the strains transmitted: locusTerm() gives log f_l(m) for one locus, strainCountTable() gives log P(m | n), and
// strainSummedLogLikelihood() sums loci and marginalises m. With the Latent parent summed out (see
// doc/marginal-latent-parent.md), f_l depends on the Latent parent's genotype only through its overlap with C_l:
// fTable() tabulates f for every overlap, contract() folds the Source's overlap distribution in to give the locus
// table log g_l(k, m), and latentSummedLogLikelihood() marginalises the Latent parent's COI k and m.
//
// The arithmetic, including operation order and the zero-probability cutoffs, is the code it was moved from, so
// results are bit-identical.
//
// The transmission leak ε (0 by default): the parent pool at a locus becomes (1 - ε) p + ε π, π the locus' Allele
// frequencies, so a child allele no parent carries is improbable rather than impossible. With ε = 0 the kernel takes
// the code paths above unchanged; with ε > 0 f is computed from the leaky pool by a generating function with no
// cancellation (leakyLocusTerms()). See section 6 of doc/marginal-latent-parent.md.
//

#ifndef TRANSMISSION_NETWORKS_APP_TRANSMISSIONKERNEL_H
#define TRANSMISSION_NETWORKS_APP_TRANSMISSIONKERNEL_H

#include "core/datatypes/Genotype.h"
#include "core/utils/ProbAnyMissing.h"
#include "core/utils/fastmath.h"
#include "core/utils/numerics.h"
#include "model/transmission_process/source_transmission_process/SourceOverlap.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace transmission_nets::model::transmission_process {

    template<int MAX_PARENTS, int MAX_STRAINS, int MAX_COI>
    class TransmissionKernel {
    public:
        using Genotype = core::datatypes::Genotype;

    public:
        // Throws std::invalid_argument unless 0 <= `leak` < 1.
        explicit TransmissionKernel(const double leak = 0) : leak_(checkedLeak(leak)) {}

        // The transmission leak ε.
        [[nodiscard]] double leak() const { return leak_; }

        using Likelihood  = double;
        using Probability = double;
        static constexpr Likelihood NEG_INF = -std::numeric_limits<Likelihood>::infinity();
        // Cells of a latent locus table: k * MAX_STRAINS + (m - 1); k = 0 is unused.
        static constexpr int NUM_CELLS = (MAX_COI + 1) * MAX_STRAINS;

        // How a locus' parent allele pool is built: from the observed parents, from them plus a latent parent, or from
        // a latent parent alone.
        enum class Variant : unsigned char {
            ObservedParents,
            ObservedAndLatentParents,
            LatentParentOnly
        };

        using StrainLogLikelihoods = std::array<Likelihood, MAX_STRAINS>;

        // Contribution of one locus to each strain count. `impossible` means the whole parent set has zero likelihood.
        struct LocusTerm {
            bool impossible = false;
            StrainLogLikelihoods terms{};
        };

        // Locus terms summed over a child's loci, kept so that one locus' term can be swapped for another: -inf terms
        // and impossible loci are counted rather than added, so removing them is exact.
        struct StrainSums {
            StrainLogLikelihoods finite{};
            std::array<std::uint16_t, MAX_STRAINS> infinite{};
            std::uint16_t impossible = 0;
            // Set when summing stopped at this locus, the first that makes the parent set impossible. The sums then
            // cover only the loci before it, so no term can be swapped in or out; the parent set stays impossible
            // until this locus' term changes.
            const void* stoppedAt = nullptr;

            void add(const LocusTerm& term) {
                if (term.impossible) {
                    ++impossible;
                    return;
                }
                for (int idx = 0; idx < MAX_STRAINS; ++idx) {
                    if (term.terms[idx] == NEG_INF) {
                        ++infinite[idx];
                    } else {
                        finite[idx] += term.terms[idx];
                    }
                }
            }
            void remove(const LocusTerm& term) {
                if (term.impossible) {
                    --impossible;
                    return;
                }
                for (int idx = 0; idx < MAX_STRAINS; ++idx) {
                    if (term.terms[idx] == NEG_INF) {
                        --infinite[idx];
                    } else {
                        finite[idx] -= term.terms[idx];
                    }
                }
            }
        };

        // log P(m | n) for n = 1..MAX_PARENTS + 1 parents and m = 1..MAX_STRAINS strains: entry (n - 1) * MAX_STRAINS
        // + (m - 1). The strains transmitted from one parent are zero-truncated Poisson with rate lambda; from several,
        // their sum, with fewer strains than parents impossible.
        using StrainCountTable = std::array<Likelihood, MAX_STRAINS * (MAX_PARENTS + 1)>;

        static StrainCountTable strainCountTable(const double lambda) {
            StrainCountTable table;
            table.fill(NEG_INF);
            for (int kk = 0; kk < MAX_PARENTS + 1; ++kk) {
                const int num_parents   = kk + 1;
                Probability denominator = 0.0;
                for (int jj = kk; jj < MAX_STRAINS; ++jj) {
                    const int num_strains = jj + 1;
                    double correction     = 0;
                    for (int ii = 0; ii < num_parents; ++ii) {
                        correction += std::pow(-1, ii) * std::pow(num_parents - ii, num_strains) * boost::math::binomial_coefficient<double>(num_parents, ii);
                    }
                    table[kk * MAX_STRAINS + jj] = num_strains * std::log(lambda) - (num_parents * std::log(std::exp(lambda) - 1)) - std::log(boost::math::factorial<double>(num_strains)) + std::log(correction);
                    denominator += std::exp(table[kk * MAX_STRAINS + jj]);
                }
                // Normalize
                denominator = std::log(denominator);
                for (int jj = 0; jj < MAX_STRAINS; ++jj) {
                    table[kk * MAX_STRAINS + jj] -= denominator;
                }
            }
            return table;
        }

        static Likelihood probNumStrains(const StrainCountTable& table, const int num_strains, const int num_parents) {
            return table[(num_parents - 1) * MAX_STRAINS + (num_strains - 1)];
        }

        // log f_l(m) for one locus from the genotypes at it: the child's, the observed parents' in parent-set order,
        // and the latent parent's (null for the observed-parents variant). `frequencies` are the locus' Allele
        // frequencies, by allele index; read only with a leak, and required then.
        LocusTerm locusTerm(const Variant variant, const Genotype& childGenotype, const std::span<const Genotype* const> parents, const Genotype* latentParent, const std::vector<double>* frequencies = nullptr) {
            LocusTerm out{};
            const size_t numParents = parents.size() + (latentParent != nullptr ? 1 : 0);

            // Only the child's alleles' frequencies are used, so only they are accumulated, in index order; each receives
            // the same additions, in the same parent order, as when every allele was accumulated.
            childAlleles_.clear();
            for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
                if (childGenotype.allele(i)) {
                    childAlleles_.push_back(i);
                }
            }
            parentPopFreqs_.assign(childAlleles_.size(), 0.0);

            double zeroProbCutoff = 1e-6;
            switch (variant) {
                case Variant::ObservedParents: {
                    zeroProbCutoff = 1e-10;
                    for (const Genotype* parent : parents) {
                        const auto& parentGenotype    = *parent;
                        const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                        for (size_t k = 0; k < childAlleles_.size(); ++k) {
                            if (parentGenotype.allele(childAlleles_[k])) {
                                parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / (totalAllelesPresent * numParents));
                            }
                        }
                    }
                    break;
                }
                case Variant::ObservedAndLatentParents: {
                    for (const Genotype* parent : parents) {
                        const auto& parentGenotype    = *parent;
                        const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                        for (size_t k = 0; k < childAlleles_.size(); ++k) {
                            if (parentGenotype.allele(childAlleles_[k])) {
                                parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / totalAllelesPresent / static_cast<Probability>(numParents));
                            }
                        }
                    }
                    const auto& latentParentGenotype = *latentParent;
                    const int totalAllelesPresent    = latentParentGenotype.totalPositiveCount();
                    if (Genotype::truePositiveCount(latentParentGenotype, childGenotype) == 0) {
                        out.impossible = true;
                        return out;
                    }
                    for (size_t k = 0; k < childAlleles_.size(); ++k) {
                        parentPopFreqs_[k] += (static_cast<Probability>(latentParentGenotype.allele(childAlleles_[k])) / totalAllelesPresent / static_cast<Probability>(numParents));
                    }
                    break;
                }
                case Variant::LatentParentOnly: {
                    const auto& parentGenotype    = *latentParent;
                    const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                    for (size_t k = 0; k < childAlleles_.size(); ++k) {
                        if (parentGenotype.allele(childAlleles_[k])) {
                            parentPopFreqs_[k] += (static_cast<Probability>(parentGenotype.allele(childAlleles_[k])) / totalAllelesPresent);
                        }
                    }
                    break;
                }
            }

            if (leak_ > 0) {
                // The leaky pool on the child's alleles: (1 - ε) p_j + ε π_j.
                if (frequencies == nullptr) {
                    throw std::logic_error("a leaky locus term needs the locus' Allele frequencies");
                }
                for (size_t k = 0; k < childAlleles_.size(); ++k) {
                    parentPopFreqs_[k] = (1.0 - leak_) * parentPopFreqs_[k] + leak_ * (*frequencies)[childAlleles_[k]];
                }
                leakyLocusTerms(out);
                return out;
            }

            Probability constrainedSetProb = 0.0;
            prVec_.clear();
            bool zeroProbEvent = false;
            for (size_t k = 0; k < childAlleles_.size(); ++k) {
                prVec_.push_back(parentPopFreqs_[k]);
                constrainedSetProb += parentPopFreqs_[k];
                zeroProbEvent = zeroProbEvent || std::abs(parentPopFreqs_[k]) < zeroProbCutoff;
            }
            if (prVec_.empty() || zeroProbEvent) {
                out.impossible = true;
                return out;
            }
            for (Likelihood& af : prVec_) {
                af /= constrainedSetProb;
            }
            const Likelihood logConstrainedSetProb   = std::log(constrainedSetProb);
            const std::vector<Likelihood>& pamVec = probAnyMissing_.vectorized(prVec_, MAX_STRAINS);
            for (int numStrains = 1; numStrains <= MAX_STRAINS; ++numStrains) {
                const int idx   = numStrains - 1;
                out.terms[idx] = pamVec[idx] >= 1.0 ? NEG_INF : std::log(1.0 - pamVec[idx]) + logConstrainedSetProb * numStrains;
            }
            return out;
        }

        // log f(m) from the pool on the child's alleles in parentPopFreqs_ (the pool's other alleles enter only through
        // these, as they sum to 1 with them): the probability that m draws from the pool show exactly the child's
        // alleles is m! [t^m] prod_{j in C} (e^{p_j t} - 1). Every coefficient is non-negative, so there is no
        // cancellation however small a leaked p_j is.
        void leakyLocusTerms(LocusTerm& out) const {
            const size_t c = parentPopFreqs_.size();
            if (c == 0 or c > static_cast<size_t>(MAX_STRAINS)) {
                out.impossible = true;
                return;
            }
            std::array<double, MAX_STRAINS + 1> series{};
            std::array<double, MAX_STRAINS + 1> next{};
            series[0] = 1.0;
            std::array<double, MAX_STRAINS + 1> drawn{};
            for (size_t k = 0; k < c; ++k) {
                // e^{p t} - 1 to degree MAX_STRAINS.
                const double p = parentPopFreqs_[k];
                double power   = 1.0;
                for (int b = 1; b <= MAX_STRAINS; ++b) {
                    power *= p;
                    drawn[b] = power * INVERSE_FACTORIAL[b];
                }
                // Each allele so far needs a draw, so series[a] is 0 for a < k.
                next.fill(0.0);
                for (int a = static_cast<int>(k); a < MAX_STRAINS; ++a) {
                    if (series[a] == 0) {
                        continue;
                    }
                    for (int b = 1; a + b <= MAX_STRAINS; ++b) {
                        next[a + b] += series[a] * drawn[b];
                    }
                }
                series = next;
            }
            bool any = false;
            for (int m = 1; m <= MAX_STRAINS; ++m) {
                const double f   = series[m] * FACTORIAL[m];
                out.terms[m - 1] = f > 0 ? std::log(f) : NEG_INF;
                any              = any or f > 0;
            }
            out.impossible = not any;
        }

        // A parent set's log-likelihood marginalised over strain counts, before the parent-set-size term. numParents
        // counts the latent parent.
        static Likelihood strainSummedLogLikelihood(const StrainSums& sums, const size_t numParents, const StrainCountTable& table) {
            if (sums.impossible) {
                return NEG_INF;
            }
            std::array<Likelihood, MAX_STRAINS> logLikelihoods{};
            for (int idx = 0; idx < MAX_STRAINS; ++idx) {
                logLikelihoods[idx] = sums.infinite[idx] > 0 ? NEG_INF : sums.finite[idx];
            }
            // Every parent transmits at least one strain, so fewer strains than parents is impossible. These entries
            // used to enter the log-sum-exp without a strain-count probability, inflating multi-parent sets.
            for (size_t numStrains = 1; numStrains < numParents; ++numStrains) {
                logLikelihoods[numStrains - 1] = NEG_INF;
            }
            Likelihood maxLl = NEG_INF;
            for (size_t numStrains = numParents; numStrains <= static_cast<size_t>(MAX_STRAINS); ++numStrains) {
                const size_t idx = numStrains - 1;
                // Add the probability of the number of strains
                if (logLikelihoods[idx] == NEG_INF) {
                    continue;
                }
                logLikelihoods[idx] += probNumStrains(table, static_cast<int>(numStrains), static_cast<int>(numParents));
                maxLl = std::max(maxLl, logLikelihoods[idx]);
            }
            return core::utils::logSumExpKnownMax(logLikelihoods.begin(), logLikelihoods.end(), maxLl);
        }

        // log g(k, m) at one locus with the latent parent summed out. Where g(k, m) > 0 is a rectangle: k >= kMin (the
        // fewest alleles a latent parent needs for the parent set to be possible) and m >= mMin (the child's alleles).
        // Cells outside it hold 0 and are excluded by the bounds, so tables add as plain arrays. `impossible` when no
        // cell is positive.
        struct LocusTable {
            // Leaves `terms` uninitialised: tables are built in place in a memo, which would otherwise zero them.
            LocusTable() noexcept {}
            bool impossible   = false;
            std::uint8_t kMin = 0;
            std::uint8_t mMin = 0;
            std::array<Likelihood, NUM_CELLS> terms;
        };

        // Locus tables summed over a child's loci, kept so that one locus' table can be swapped for another. The sum's
        // support is the intersection of the tables' rectangles, so the bounds are kept as counts: removing a table is
        // exact.
        struct LatentSums {
            std::array<Likelihood, NUM_CELLS> finite{};
            std::array<std::uint16_t, MAX_COI + 1> kMinCount{};
            std::array<std::uint16_t, MAX_STRAINS + 1> mMinCount{};
            std::uint16_t impossible = 0;
            // As StrainSums::stoppedAt.
            const void* stoppedAt = nullptr;

            void add(const LocusTable& table) {
                if (table.impossible) {
                    ++impossible;
                    return;
                }
                const int first = table.kMin * MAX_STRAINS + (table.mMin - 1);
                for (int i = first; i < NUM_CELLS; ++i) {
                    finite[i] += table.terms[i];
                }
                ++kMinCount[table.kMin];
                ++mMinCount[table.mMin];
            }
            void remove(const LocusTable& table) {
                if (table.impossible) {
                    --impossible;
                    return;
                }
                const int first = table.kMin * MAX_STRAINS + (table.mMin - 1);
                for (int i = first; i < NUM_CELLS; ++i) {
                    finite[i] -= table.terms[i];
                }
                --kMinCount[table.kMin];
                --mMinCount[table.mMin];
            }
            // The support's lower bounds.
            [[nodiscard]] int kLow() const {
                for (int k = MAX_COI; k > 0; --k) {
                    if (kMinCount[k] > 0) {
                        return k;
                    }
                }
                return 1;
            }
            [[nodiscard]] int mLow() const {
                for (int m = MAX_STRAINS; m > 0; --m) {
                    if (mMinCount[m] > 0) {
                        return m;
                    }
                }
                return 1;
            }
        };

        // f(m; U, w) for the nonempty U ⊆ C (a mask over C's alleles in index order) and w = |G \ C| up to maxW for
        // which the parent set is possible, as probabilities; only those are kept. With the bounds of where g can be
        // positive: k >= rMin, the fewest alleles of a possible latent parent, and m >= mMin.
        struct FEntry {
            std::uint32_t mask = 0;
            std::uint16_t w    = 0;
            std::uint16_t r    = 0;// |U| + w
            std::array<double, MAX_STRAINS> f{};
        };
        struct FTable {
            unsigned int maxW = 0;
            std::vector<FEntry> entries{};
            int rMin = MAX_COI + 1;
            int mMin = MAX_STRAINS + 1;
        };

        // The F table for a child genotype and the observed parents' genotypes (none for the latent parent alone).
        // `frequencies` as for locusTerm(): with a leak the table depends on them.
        FTable fTable(const Genotype& child, const std::span<const Genotype* const> parents, const std::vector<double>* frequencies = nullptr) {
            if (leak_ > 0 and frequencies == nullptr) {
                throw std::logic_error("a leaky F table needs the locus' Allele frequencies");
            }
            const auto childAlleles = allelesOf(child);
            const unsigned int A    = child.totalAlleles();
            const auto c            = static_cast<unsigned int>(childAlleles.size());
            FTable table{};
            if (c == 0 or c > static_cast<unsigned int>(MAX_STRAINS)) {
                // No alleles, or more than the most strains transmitted can carry: impossible.
                return table;
            }
            // A child allele in no observed parent has frequency 0 unless the latent parent carries it, so masks
            // without all of them are impossible and skipped without computing their terms. With a leak no allele is
            // required.
            unsigned int required = 0;
            for (unsigned int i = 0; i < c and leak_ == 0; ++i) {
                if (std::ranges::none_of(parents, [&](const Genotype* parent) { return parent->allele(childAlleles[i]); })) {
                    required |= 1u << i;
                }
            }
            table.maxW                  = std::min<unsigned int>(A - c, MAX_COI);
            const unsigned int numMasks = 1u << c;
            const Variant variant       = parents.empty() ? Variant::LatentParentOnly : Variant::ObservedAndLatentParents;
            // The representative latent genotypes (as representative() in the tests' MarginalLatentParentReference.h
            // builds them), built incrementally: U's alleles, then one more allele outside C, in index order, for each w.
            std::vector<unsigned int> outside;
            for (unsigned int a = 0; a < A and outside.size() < table.maxW; ++a) {
                if (not child.allele(a)) {
                    outside.push_back(a);
                }
            }
            Genotype latent = child;
            // With a leak, the latent parent alone may share no allele with the child (U = ∅): the child's alleles
            // then all come from the leak. With observed parents U = ∅ stays impossible (the latent parent must share
            // an allele with the child), as does any U = ∅ without a leak.
            const unsigned int firstMask = leak_ > 0 and parents.empty() ? 0 : 1;
            for (unsigned int mask = firstMask; mask < numMasks; ++mask) {
                if ((mask & required) != required) {
                    continue;
                }
                latent.reset();
                for (unsigned int i = 0; i < c; ++i) {
                    if (mask >> i & 1u) {
                        latent.set(childAlleles[i]);
                    }
                }
                for (unsigned int w = 0; w <= table.maxW; ++w) {
                    if (std::popcount(mask) + w > static_cast<unsigned int>(MAX_COI)) {
                        break;
                    }
                    if (w > 0) {
                        latent.set(outside[w - 1]);
                    }
                    if (mask == 0 and w == 0) {
                        continue;// a latent parent carries at least one allele
                    }
                    const auto term = locusTerm(variant, child, parents, &latent, frequencies);
                    if (term.impossible) {
                        continue;
                    }
                    FEntry entry{.mask = mask, .w = static_cast<std::uint16_t>(w), .r = static_cast<std::uint16_t>(std::popcount(mask) + w)};
                    bool any = false;
                    for (int m = 0; m < MAX_STRAINS; ++m) {
                        entry.f[m] = std::exp(term.terms[m]);
                        any        = any or entry.f[m] > 0;
                    }
                    if (any) {
                        table.rMin = std::min(table.rMin, static_cast<int>(entry.r));
                        for (int m = 0; m < MAX_STRAINS; ++m) {
                            if (entry.f[m] > 0) {
                                table.mMin = std::min(table.mMin, m + 1);
                                break;
                            }
                        }
                        table.entries.push_back(entry);
                    }
                }
            }
            table.entries.shrink_to_fit();
            return table;
        }

        // The locus table log g(k, m): the F table contracted with the Source's overlap distribution at the locus.
        static void contract(const FTable& f, const marginal::SourceOverlap& overlap, LocusTable& out) {
            if (f.entries.empty()) {
                out.impossible = true;
                return;
            }
            // Every frequency is positive, so P(U, w | k) > 0 for k >= |U| + w, and g > 0 on k >= rMin, m >= mMin.
            const int kMin = f.rMin;
            const int mMin = f.mMin;
            std::array<double, NUM_CELLS> g;
            std::fill(g.begin() + kMin * MAX_STRAINS, g.end(), 0.0);
            for (const auto& entry : f.entries) {
                const auto& fm   = entry.f;
                const double* ps = overlap.row(entry.mask, entry.w);
                for (int k = static_cast<int>(entry.r); k <= MAX_COI; ++k) {
                    const double p = ps[k];
                    double* row    = g.data() + k * MAX_STRAINS;
                    for (int m = 0; m < MAX_STRAINS; ++m) {
                        row[m] += p * fm[m];
                    }
                }
            }
            out.impossible = false;
            out.kMin       = static_cast<std::uint8_t>(kMin);
            out.mMin       = static_cast<std::uint8_t>(mMin);
            std::fill(out.terms.begin(), out.terms.end(), 0.0);
            // Inside the rectangle every cell is positive in exact arithmetic; one that underflowed below the smallest
            // normal double counts as that, which is as good as zero.
            for (int k = kMin; k <= MAX_COI; ++k) {
                const double* in = g.data() + k * MAX_STRAINS + mMin - 1;
                double* to       = out.terms.data() + k * MAX_STRAINS + mMin - 1;
                for (int j = 0; j <= MAX_STRAINS - mMin; ++j) {
                    to[j] = core::utils::logPositive(std::max(in[j], std::numeric_limits<double>::min()));
                }
            }
        }

        // The parent set's log-likelihood with the latent parent summed out, before the parent-set-size term: the sum
        // over its COI k, with log probabilities coiLogProbs[k], and the strains transmitted m. numParents counts the
        // latent parent.
        static Likelihood latentSummedLogLikelihood(const LatentSums& sums, const size_t numParents, const std::array<double, MAX_COI + 1>& coiLogProbs, const StrainCountTable& table) {
            if (sums.impossible) {
                return NEG_INF;
            }
            const int mFirst = std::max(sums.mLow(), static_cast<int>(numParents));
            const int kFirst = sums.kLow();
            std::array<Likelihood, MAX_STRAINS + 1> strainLogProbs{};
            for (int m = mFirst; m <= MAX_STRAINS; ++m) {
                strainLogProbs[m] = probNumStrains(table, m, static_cast<int>(numParents));
            }
            // The support is the rectangle k >= kFirst, m >= mFirst; the COI distribution is positive for every k >= 1.
            const int width = MAX_STRAINS - mFirst + 1;
            std::array<Likelihood, NUM_CELLS> terms;
            int count          = 0;
            Likelihood maxTerm = NEG_INF;
            for (int k = kFirst; k <= MAX_COI; ++k) {
                const Likelihood* row = sums.finite.data() + k * MAX_STRAINS + mFirst - 1;
                const Likelihood* tx  = strainLogProbs.data() + mFirst;
                const Likelihood coi  = coiLogProbs[k];
                Likelihood* out       = terms.data() + count;
                for (int j = 0; j < width; ++j) {
                    out[j]  = coi + tx[j] + row[j];
                    maxTerm = std::max(maxTerm, out[j]);
                }
                count += width;
            }
            if (maxTerm == NEG_INF) {
                return NEG_INF;
            }
            // Cells more than NEGLIGIBLE below the largest add less than exp(-NEGLIGIBLE) each; usually most are.
            constexpr Likelihood NEGLIGIBLE = 40;
            std::array<Likelihood, NUM_CELLS> near;
            int numNear = 0;
            for (int i = 0; i < count; ++i) {
                const Likelihood d = terms[i] - maxTerm;
                near[numNear]      = d;
                numNear += d > -NEGLIGIBLE ? 1 : 0;
            }
            Likelihood total = 0;
            for (int i = 0; i < numNear; ++i) {
                total += core::utils::expNonPositive(near[i]);
            }
            return maxTerm + std::log(total);
        }

        // The indices of a genotype's alleles.
        static std::vector<unsigned int> allelesOf(const Genotype& g) {
            std::vector<unsigned int> out;
            for (unsigned int i = 0; i < g.totalAlleles(); ++i) {
                if (g.allele(i)) {
                    out.push_back(i);
                }
            }
            return out;
        }

    private:
        static double checkedLeak(const double leak) {
            if (not(leak >= 0 and leak < 1)) {
                throw std::invalid_argument("the transmission leak must be at least 0 and below 1");
            }
            return leak;
        }

        static constexpr std::array<double, MAX_STRAINS + 1> FACTORIAL = [] {
            std::array<double, MAX_STRAINS + 1> out{};
            out[0] = 1.0;
            for (int b = 1; b <= MAX_STRAINS; ++b) {
                out[b] = out[b - 1] * b;
            }
            return out;
        }();
        static constexpr std::array<double, MAX_STRAINS + 1> INVERSE_FACTORIAL = [] {
            std::array<double, MAX_STRAINS + 1> out{};
            for (int b = 0; b <= MAX_STRAINS; ++b) {
                out[b] = 1.0 / FACTORIAL[b];
            }
            return out;
        }();

        double leak_ = 0;
        core::utils::probAnyMissingFunctor probAnyMissing_{};
        // Scratch buffers reused across locus evaluations.
        std::vector<size_t> childAlleles_{};
        std::vector<Probability> parentPopFreqs_{};
        std::vector<Likelihood> prVec_{};
    };

}// namespace transmission_nets::model::transmission_process

#endif//TRANSMISSION_NETWORKS_APP_TRANSMISSIONKERNEL_H
