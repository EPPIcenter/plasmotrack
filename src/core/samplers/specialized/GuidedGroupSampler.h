//
// Moves a family of infections that can infect one another to new infection times, the family's genotypes following,
// in one Metropolis-Hastings step: GuidedPairSampler's swap generalised from two infections to a few.
//

#ifndef TRANSMISSION_NETWORKS_APP_GUIDEDGROUPSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_GUIDEDGROUPSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/LogWeightedDraw.h"
#include "core/samplers/Proposal.h"
#include "core/samplers/specialized/GuidedPairSampler.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    /*
     * One step picks a group from a fixed list and proposes:
     * 1. new infection times for `movers` of its members (distinct, in a uniform order), each drawn after the previous
     *    one moved: with probability 1/2 uniformly over its feasible window, otherwise within `jitter` of another
     *    member's time (picked uniformly). The other members' times stay, so the density both ways is known.
     * 2. the group's genotypes at its units, each allele observed in any member at a locus where the members' observed
     *    genotypes differ. A deterministic pilot first makes the genotypes consistent with the best transmission tree
     *    inside the group (pilot()); a walk then draws every unit, in a fixed order, among the 2^k ways the members may
     *    hold it, in proportion to the target.
     * The reverse move puts the times back in reverse order, runs the pilot from there, and takes the walk's probability
     * of drawing the current genotypes; a reverse walk that ends anywhere else cannot reach the current state and
     * rejects. The pilot is a deterministic function of the state it starts from, so the time densities and the walk
     * probabilities both ways are the whole proposal density.
     */
    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    class GuidedGroupSampler : public AbstractSampler {
    public:
        using Member = typename GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::Member;
        using Group  = std::vector<Member>;

        GuidedGroupSampler(std::vector<Group> groups, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double minDuration, unsigned int maxCOI, double jitter = 10.0,
                           size_t movers = 2, size_t moves = 1);

        // `moves` steps.
        void update() noexcept override {
            for (size_t m = 0; m < moves_; ++m) step();
        }

        [[nodiscard]] size_t groups() const noexcept { return prepared_.size(); }

    private:
        using Genotype  = std::remove_cvref_t<decltype(std::declval<InfectionEventImpl&>().latentGenotype(std::declval<std::shared_ptr<LocusImpl>>())->value())>;
        using Option    = std::vector<Genotype>;// every member's genotype at a unit's locus
        using Genotypes = std::vector<Option>;

        struct Unit {
            std::shared_ptr<LocusImpl> locus;
            unsigned int allele;
        };
        struct Prepared {
            Group group;
            std::vector<Unit> units;
        };

        void optionsAt(const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept;
        void set(const Prepared& p, const Unit& u, const Option& o) noexcept;
        void pilot(const Prepared& p) noexcept;
        double walk(const Prepared& p, const Genotypes* towards) noexcept;
        [[nodiscard]] Genotypes snapshot(const Prepared& p) const;
        void restore(const Prepared& p, const std::vector<double>& durations, const Genotypes& genotypes) noexcept;

        std::vector<Prepared> prepared_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double minDuration_;
        unsigned int maxCOI_;
        double jitter_ = 10.0;
        size_t movers_ = 2;
        size_t moves_  = 1;
        SamplerCounts innerCounts_{};

        void step() noexcept;

        // Relocates `who` in turn to `to`; the log density of the forward draws and of the reverse ones (undoing in
        // reverse order), or -inf when a time is out of bounds or the target infinite (the state then partly moved).
        std::pair<double, double> relocateAll(const Prepared& p, const std::vector<size_t>& who, const std::vector<double>& to) noexcept;
        double drawTime(const Prepared& p, size_t i) noexcept;

        // The density of infection time t for member i under the relocation proposal, the others' times as they are.
        [[nodiscard]] double timeDensity(const Prepared& p, size_t i, double t) const noexcept;
        bool setTime(const Prepared& p, size_t i, double t) noexcept;
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::GuidedGroupSampler(std::vector<Group> groups, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, const double minDuration,
                                                                                      const unsigned int maxCOI, const double jitter, const size_t movers, const size_t moves)
        : target_(std::move(target)), rng_(std::move(rng)), minDuration_(minDuration), maxCOI_(maxCOI), jitter_(jitter), movers_(movers), moves_(moves) {
        for (auto& group : groups) {
            Prepared p{.group = std::move(group)};
            for (const auto& locus : p.group.front().infection->loci()) {
                std::set<unsigned int> alleles;
                std::set<std::string> seen;
                for (const auto& m : p.group) {
                    if (m.infection->observedGenotype().contains(locus)) {
                        const auto observed = m.infection->observedGenotype().at(locus)->value();
                        seen.insert(observed.serialize());
                        for (unsigned int allele = 0; allele < observed.totalAlleles(); ++allele) {
                            if (observed.allele(allele)) alleles.insert(allele);
                        }
                    }
                }
                // Loci where every member is observed alike have nothing to rearrange.
                if (seen.size() <= 1) continue;
                for (const auto allele : alleles) p.units.push_back({locus, allele});
            }
            prepared_.push_back(std::move(p));
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::optionsAt(const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept {
        const size_t k = p.group.size();
        Option current(k);
        for (size_t i = 0; i < k; ++i) current[i] = p.group[i].infection->latentGenotype(u.locus)->value();
        options.clear();
        values.clear();
        for (unsigned int bits = 0; bits < (1u << k); ++bits) {
            Option o = current;
            bool ok  = true;
            for (size_t i = 0; i < k and ok; ++i) {
                o[i].set(u.allele, bits >> i & 1u);
                ok = o[i].totalPositiveCount() > 0 and o[i].totalPositiveCount() <= maxCOI_;
            }
            if (not ok) continue;
            Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
            for (size_t i = 0; i < k; ++i) {
                if (o[i] != current[i]) move.change(*p.group[i].infection->latentGenotype(u.locus), o[i]);
            }
            values.push_back(move.value());
            move.revert();
            options.push_back(std::move(o));
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::set(const Prepared& p, const Unit& u, const Option& o) noexcept {
        Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
        bool any = false;
        for (size_t i = 0; i < p.group.size(); ++i) {
            auto& g = *p.group[i].infection->latentGenotype(u.locus);
            if (o[i] != g.value()) {
                move.change(g, o[i]);
                any = true;
            }
        }
        if (any) {
            move.accept();
        } else {
            move.revert();
        }
    }

    // Genotypes made consistent with a hypothesised tree inside the group. Members in their current time order; each but
    // the first gets a parent among the earlier members or none (outside the group). For each assignment, at every
    // unit's locus: a child's observed alleles its parent lacks are called false positives when at most one is and some
    // remain, otherwise the parent gains them; then hidden alleles (in the latent genotype, not observed) that no child
    // in the group carries are dropped. The best assignment's genotypes are kept, if better than the current ones.
    // Deterministic.
    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::pilot(const Prepared& p) noexcept {
        const size_t k = p.group.size();
        std::vector<size_t> byTime(k);
        std::iota(byTime.begin(), byTime.end(), 0);
        std::ranges::sort(byTime, [&](size_t a, size_t b) { return p.group[a].infection->infectionTime() < p.group[b].infection->infectionTime(); });
        std::vector<std::shared_ptr<LocusImpl>> loci;
        for (const auto& u : p.units) {
            if (std::ranges::find(loci, u.locus) == loci.end()) loci.push_back(u.locus);
        }
        // Assignment: parent position (in byTime) for positions 1..k-1, k meaning none.
        std::vector<size_t> assign(k, k);
        double bestValue = target_->value();
        std::vector<std::vector<Genotype>> bestGenotypes;
        auto build = [&](const std::vector<size_t>& a) {
            // Genotypes per member per locus under assignment a, from the current ones.
            std::vector<std::vector<Genotype>> g(k);
            for (size_t i = 0; i < k; ++i) {
                for (const auto& locus : loci) g[i].push_back(p.group[i].infection->latentGenotype(locus)->value());
            }
            for (size_t l = 0; l < loci.size(); ++l) {
                for (size_t pos = 1; pos < k; ++pos) {
                    if (a[pos] == k) continue;
                    const size_t c = byTime[pos], q = byTime[a[pos]];
                    const auto& inf = *p.group[c].infection;
                    if (not inf.observedGenotype().contains(loci[l])) continue;
                    const auto obs = inf.observedGenotype().at(loci[l])->value();
                    if (obs.totalPositiveCount() == 0) continue;
                    const auto lacking = Genotype::minus(obs, g[q][l]);// observed in the child, absent in the parent
                    const auto kept    = Genotype::shared(obs, g[q][l]);
                    if (lacking.totalPositiveCount() == 0) {
                        g[c][l] = Genotype::any(g[c][l], obs);
                    } else if (lacking.totalPositiveCount() <= 1 and kept.totalPositiveCount() > 0) {
                        g[c][l] = kept;
                    } else {
                        g[q][l] = Genotype::any(g[q][l], obs);
                        g[c][l] = Genotype::any(Genotype::shared(g[c][l], g[q][l]), obs);
                    }
                }
                // Drop hidden alleles no group child carries.
                for (size_t i = 0; i < k; ++i) {
                    const auto& inf = *p.group[i].infection;
                    Genotype observed = inf.observedGenotype().contains(loci[l]) ? inf.observedGenotype().at(loci[l])->value() : g[i][l];
                    if (observed.totalPositiveCount() == 0) continue;
                    auto hidden = Genotype::minus(g[i][l], observed);
                    if (hidden.totalPositiveCount() == 0) continue;
                    Genotype needed = observed;
                    for (size_t pos = 1; pos < k; ++pos) {
                        if (a[pos] != k and byTime[a[pos]] == i) needed = Genotype::any(needed, Genotype::shared(g[i][l], g[byTime[pos]][l]));
                    }
                    g[i][l] = Genotype::shared(g[i][l], Genotype::any(observed, needed));
                    if (g[i][l].totalPositiveCount() == 0) g[i][l] = observed;
                }
            }
            return g;
        };
        // Puts genotypes g (per member, per locus) into a Proposal, leaving out any outside 1 to maxCOI Alleles.
        auto change = [&](Proposal<T, Engine>& move, const std::vector<std::vector<Genotype>>& g) {
            for (size_t i = 0; i < k; ++i) {
                for (size_t l = 0; l < loci.size(); ++l) {
                    auto& param = *p.group[i].infection->latentGenotype(loci[l]);
                    if (g[i][l].totalPositiveCount() > 0 and g[i][l].totalPositiveCount() <= maxCOI_ and param.value() != g[i][l]) move.change(param, g[i][l]);
                }
            }
        };
        while (true) {
            const auto g = build(assign);
            Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
            change(move, g);
            const double v = move.value();
            move.revert();
            if (v > bestValue) {
                bestValue     = v;
                bestGenotypes = g;
            }
            // Next assignment: position pos ranges over 0..pos-1 and k.
            size_t pos = 1;
            for (; pos < k; ++pos) {
                size_t& x = assign[pos];
                x = x == k ? 0 : x + 1;
                if (x == pos) x = k;// past the earlier positions: wrap to none and carry
                if (x != k) break;
            }
            if (pos == k) break;
        }
        if (not bestGenotypes.empty()) {
            Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
            change(move, bestGenotypes);
            move.value();
            move.accept();
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    double GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::walk(const Prepared& p, const Genotypes* towards) noexcept {
        std::vector<Option> options;
        std::vector<double> values;
        double logQ = 0;
        for (size_t u = 0; u < p.units.size(); ++u) {
            const auto allele = p.units[u].allele;
            optionsAt(p, p.units[u], options, values);
            if (options.empty()) return -std::numeric_limits<double>::infinity();
            const auto weights = logWeights(values);
            if (not std::isfinite(weights.top)) return -std::numeric_limits<double>::infinity();
            size_t pick = 0;
            if (towards != nullptr) {
                const auto& want = (*towards)[u];
                pick = static_cast<size_t>(std::ranges::find_if(options, [&](const Option& o) {
                                               for (size_t i = 0; i < o.size(); ++i) {
                                                   if (o[i].allele(allele) != want[i].allele(allele)) return false;
                                               }
                                               return true;
                                           }) - options.begin());
                if (pick == options.size()) return -std::numeric_limits<double>::infinity();
            } else {
                pick = drawLogWeighted(values, weights, *rng_);
            }
            logQ += values[pick] - weights.top - std::log(weights.total);
            set(p, p.units[u], options[pick]);
        }
        return logQ;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    typename GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::Genotypes GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::snapshot(const Prepared& p) const {
        Genotypes out;
        for (const auto& u : p.units) {
            Option o;
            for (const auto& m : p.group) o.push_back(m.infection->latentGenotype(u.locus)->value());
            out.push_back(std::move(o));
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::restore(const Prepared& p, const std::vector<double>& durations, const Genotypes& genotypes) noexcept {
        Proposal<T, Engine> back(*target_, *rng_, innerCounts_);
        for (size_t i = 0; i < p.group.size(); ++i) {
            if (p.group[i].infection->infectionDuration()->value() != durations[i]) back.change(*p.group[i].infection->infectionDuration(), durations[i]);
        }
        std::set<const void*> done;
        for (size_t u = 0; u < p.units.size(); ++u) {
            if (not done.insert(p.units[u].locus.get()).second) continue;
            for (size_t i = 0; i < p.group.size(); ++i) {
                auto& g = *p.group[i].infection->latentGenotype(p.units[u].locus);
                if (g.value() != genotypes[u][i]) back.change(g, genotypes[u][i]);
            }
        }
        back.value();
        back.accept();
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    double GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::timeDensity(const Prepared& p, const size_t i, const double t) const noexcept {
        const auto& m    = p.group[i];
        const double obs = m.infection->observationTime()->value();
        const double lo = obs - m.maxDuration, hi = obs - minDuration_;
        double d = (t >= lo and t <= hi) ? 0.5 / (hi - lo) : 0.0;
        const double others = static_cast<double>(p.group.size() - 1);
        for (size_t j = 0; j < p.group.size(); ++j) {
            if (j == i) continue;
            if (std::abs(t - p.group[j].infection->infectionTime()) <= jitter_) d += 0.5 / others / (2 * jitter_);
        }
        return d;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    bool GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::setTime(const Prepared& p, const size_t i, const double t) noexcept {
        const auto& m         = p.group[i];
        const double duration = m.infection->observationTime()->value() - t;
        if (duration < minDuration_ or duration > m.maxDuration) return false;
        Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
        move.change(*m.infection->infectionDuration(), duration);
        if (not std::isfinite(move.value())) {
            move.revert();
            return false;
        }
        move.accept();
        return true;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    double GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::drawTime(const Prepared& p, const size_t i) noexcept {
        boost::random::uniform_01<> u01;
        const auto& m    = p.group[i];
        const double obs = m.infection->observationTime()->value();
        const size_t k   = p.group.size();
        if (k == 1 or u01(*rng_) < 0.5) return obs - m.maxDuration + u01(*rng_) * (m.maxDuration - minDuration_);
        size_t j = boost::random::uniform_int_distribution<size_t>(0, k - 2)(*rng_);
        if (j >= i) ++j;
        return p.group[j].infection->infectionTime() + (2 * u01(*rng_) - 1) * jitter_;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    std::pair<double, double> GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::relocateAll(const Prepared& p, const std::vector<size_t>& who, const std::vector<double>& to) noexcept {
        double logF = 0, logR = 0;
        for (size_t s = 0; s < who.size(); ++s) {
            const size_t i     = who[s];
            const double from  = p.group[i].infection->infectionTime();
            const double qF = timeDensity(p, i, to[s]), qR = timeDensity(p, i, from);
            if (qF <= 0 or qR <= 0 or not setTime(p, i, to[s])) return {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
            logF += std::log(qF);
            logR += std::log(qR);
        }
        return {logF, logR};
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedGroupSampler<T, Engine, InfectionEventImpl, LocusImpl>::step() noexcept {
        if (prepared_.empty()) return;
        const auto& p  = prepared_[boost::random::uniform_int_distribution<size_t>(0, prepared_.size() - 1)(*rng_)];
        const size_t k = p.group.size();
        // An ordered sequence of distinct movers, uniform; each draws its time after the previous one moved.
        std::vector<size_t> order(k);
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), *rng_);
        const std::vector<size_t> who(order.begin(), order.begin() + static_cast<long>(std::min(movers_, k)));
        boost::random::uniform_01<> u01;
        const double before = target_->value();
        std::vector<double> durations(k);
        for (size_t a = 0; a < k; ++a) durations[a] = p.group[a].infection->infectionDuration()->value();
        const auto current = snapshot(p);
        std::vector<double> oldTimes, newTimes;
        double logTF = 0, logTR = 0;
        bool ok = true;
        for (const size_t i : who) {
            oldTimes.push_back(p.group[i].infection->infectionTime());
            const double t = drawTime(p, i);
            newTimes.push_back(t);
            const auto [f, r] = relocateAll(p, {i}, {t});
            if (not std::isfinite(f)) {
                ok = false;
                break;
            }
            logTF += f;
            logTR += r;
        }
        if (not ok) {
            restore(p, durations, current);
            countProposal(false);
            return;
        }
        pilot(p);
        const double logForward = walk(p, nullptr);
        const double after      = target_->value();
        std::vector<double> proposedDurations(k);
        for (size_t a = 0; a < k; ++a) proposedDurations[a] = p.group[a].infection->infectionDuration()->value();
        const auto proposed = snapshot(p);
        double logReverse   = -std::numeric_limits<double>::infinity();
        if (std::isfinite(logForward)) {
            // Undo the moves in reverse order: their draw densities are logTR, computed in the same surroundings.
            bool back = true;
            for (size_t s = who.size(); s-- > 0 and back;) back = setTime(p, who[s], oldTimes[s]);
            if (back) {
                pilot(p);
                logReverse = walk(p, &current);
                // A reverse walk that ends anywhere but the current state means the reverse move cannot reach it.
                if (std::isfinite(logReverse) and snapshot(p) != current) logReverse = -std::numeric_limits<double>::infinity();
            }
        }
        const double logAccept = std::isfinite(logReverse) ? after - before + logReverse - logForward + logTR - logTF : -std::numeric_limits<double>::infinity();
        const bool accepted    = std::log(u01(*rng_)) < logAccept;
        restore(p, accepted ? proposedDurations : durations, accepted ? proposed : current);
        countProposal(accepted);
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_GUIDEDGROUPSAMPLER_H
