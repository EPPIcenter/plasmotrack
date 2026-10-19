//
// Moves a pair of infections that can infect each other between orderings, the pair's genotypes following, in one
// Metropolis-Hastings step: the infection times exchanged, then the genotypes drawn where the exchange leaves them.
//

#ifndef TRANSMISSION_NETWORKS_APP_GUIDEDPAIRSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_GUIDEDPAIRSAMPLER_H

#include "core/samplers/AbstractSampler.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    /*
     * A guided swap on one pair of infections, picked uniformly from a fixed list, `moves` times per update. Changing
     * which of two infections came first can cost hundreds of log units until several alleles of both change with it.
     * The swap exchanges the pair's infection times, then sets the pair's genotypes in two steps, at the units: each
     * allele observed in either infection, which both infections' genotypes may hold or not, four options a unit.
     *
     * The pilot climbs: it scores every unit's best option once, then repeatedly scores the leading unit again and
     * sets it if it still leads, the other units at its locus scored again after, until no unit gains. The walk then
     * draws every unit in a fixed order, in proportion to the target. The reverse move is the same from the proposed
     * state: times exchanged back, the pilot, and the walk's probability of drawing the current genotypes; the swap
     * is accepted on the target's change and both walks' probabilities. The pilot is deterministic, so it needs no
     * correction of its own, and the time exchange is its own inverse.
     *
     * A genotype with no allele or more than `maxCOI` is not an option. A swap is rejected when either exchange
     * leaves a duration outside its bounds or the target at no probability; a rejected swap puts both infections
     * back as they were.
     */
    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    class GuidedPairSampler : public AbstractSampler {
    public:
        struct Member {
            std::shared_ptr<InfectionEventImpl> infection;
            double maxDuration;
        };

        struct Pair {
            Member first, second;
        };

        // Every pair of infections that can infect each other, each infection's duration bounded by `maxDuration(i)`.
        template<typename AllowedRelationships, typename MaxDuration>
        static std::vector<Pair> mutualPairs(const std::vector<std::shared_ptr<InfectionEventImpl>>& infections, AllowedRelationships& allowed, MaxDuration&& maxDuration);

        GuidedPairSampler(std::vector<Pair> pairs, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double minDuration, unsigned int maxCOI, int moves = 1);

        void update() noexcept override;

    private:
        using Genotype  = std::remove_cvref_t<decltype(std::declval<InfectionEventImpl&>().latentGenotype(std::declval<std::shared_ptr<LocusImpl>>())->value())>;
        using Option    = std::pair<Genotype, Genotype>;
        using Genotypes = std::vector<Option>;// both genotypes at each unit's locus

        struct Unit {
            std::shared_ptr<LocusImpl> locus;
            unsigned int allele;
        };

        struct Prepared {
            Pair pair;
            std::vector<Unit> units;
        };

        // Exchanges the pair's infection times if both durations stay in bounds and the target finite.
        bool exchange(const Prepared& p) noexcept;
        // The allowed options at unit u and the target at each.
        void optionsAt(const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept;
        // The best option at unit u and what it gains.
        std::pair<double, Option> bestAt(const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept;
        void set(const Prepared& p, const Unit& u, const Option& o) noexcept;
        void pilot(const Prepared& p) noexcept;
        // Draws every unit in proportion to the target, or with `towards` takes towards' alleles; the log probability.
        double walk(const Prepared& p, const Genotypes* towards) noexcept;
        [[nodiscard]] Genotypes snapshot(const Prepared& p) const;
        void restore(const Prepared& p, double durationA, double durationB, const Genotypes& genotypes) noexcept;
        void swap(const Prepared& p) noexcept;

        std::vector<Prepared> prepared_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double minDuration_;
        unsigned int maxCOI_;
        int moves_;
        SamplerCounts innerCounts_{};
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    template<typename AllowedRelationships, typename MaxDuration>
    std::vector<typename GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::Pair> GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::mutualPairs(
            const std::vector<std::shared_ptr<InfectionEventImpl>>& infections, AllowedRelationships& allowed, MaxDuration&& maxDuration) {
        std::vector<Pair> out;
        for (size_t i = 0; i < infections.size(); ++i) {
            const auto parents = allowed.permittedParents(infections[i]);
            for (size_t j = i + 1; j < infections.size(); ++j) {
                if (std::ranges::find(parents, infections[j]) == parents.end()) {
                    continue;
                }
                if (const auto theirs = allowed.permittedParents(infections[j]); std::ranges::find(theirs, infections[i]) != theirs.end()) {
                    out.push_back({{infections[i], maxDuration(i)}, {infections[j], maxDuration(j)}});
                }
            }
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::GuidedPairSampler(std::vector<Pair> pairs,
                                                                                std::shared_ptr<T> target,
                                                                                std::shared_ptr<Engine> rng,
                                                                                const double minDuration,
                                                                                const unsigned int maxCOI,
                                                                                const int moves)
        : target_(std::move(target)),
          rng_(std::move(rng)),
          minDuration_(minDuration),
          maxCOI_(maxCOI),
          moves_(moves) {
        for (auto& pair : pairs) {
            Prepared p{.pair = std::move(pair)};
            for (const auto& locus : p.pair.first.infection->loci()) {
                std::set<unsigned int> alleles;
                for (const auto* member : {&p.pair.first, &p.pair.second}) {
                    if (member->infection->observedGenotype().contains(locus)) {
                        const auto observed = member->infection->observedGenotype().at(locus)->value();
                        for (unsigned int allele = 0; allele < observed.totalAlleles(); ++allele) {
                            if (observed.allele(allele)) {
                                alleles.insert(allele);
                            }
                        }
                    }
                }
                for (const auto allele : alleles) {
                    p.units.push_back({locus, allele});
                }
            }
            if (not p.units.empty()) {
                prepared_.push_back(std::move(p));
            }
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::update() noexcept {
        if (prepared_.empty()) {
            return;
        }
        for (int m = 0; m < moves_; ++m) {
            swap(prepared_[boost::random::uniform_int_distribution<size_t>(0, prepared_.size() - 1)(*rng_)]);
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    bool GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::exchange(const Prepared& p) noexcept {
        const auto& a   = *p.pair.first.infection;
        const auto& b   = *p.pair.second.infection;
        const double da = a.observationTime()->value() - b.infectionTime(), db = b.observationTime()->value() - a.infectionTime();
        if (da < minDuration_ or da > p.pair.first.maxDuration or db < minDuration_ or db > p.pair.second.maxDuration) {
            return false;
        }
        Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
        move.change(*a.infectionDuration(), da);
        move.change(*b.infectionDuration(), db);
        if (not std::isfinite(move.value())) {
            move.revert();
            return false;
        }
        move.accept();
        return true;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::optionsAt(const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept {
        const auto& a  = *p.pair.first.infection;
        const auto& b  = *p.pair.second.infection;
        const auto ga  = a.latentGenotype(u.locus)->value(), gb = b.latentGenotype(u.locus)->value();
        auto counted   = [&](const Genotype& g) { return g.totalPositiveCount() > 0 and g.totalPositiveCount() <= maxCOI_; };
        options.clear();
        values.clear();
        for (const bool inA : {false, true}) {
            for (const bool inB : {false, true}) {
                auto na = ga, nb = gb;
                na.set(u.allele, inA);
                nb.set(u.allele, inB);
                if (not counted(na) or not counted(nb)) {
                    continue;
                }
                Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
                if (na != ga) {
                    move.change(*a.latentGenotype(u.locus), na);
                }
                if (nb != gb) {
                    move.change(*b.latentGenotype(u.locus), nb);
                }
                values.push_back(move.value());
                move.revert();
                options.emplace_back(na, nb);
            }
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    std::pair<double, typename GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::Option> GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::bestAt(
            const Prepared& p, const Unit& u, std::vector<Option>& options, std::vector<double>& values) noexcept {
        const double now = target_->value();
        optionsAt(p, u, options, values);
        const auto k = static_cast<size_t>(std::ranges::max_element(values) - values.begin());
        return {values[k] - now, options[k]};
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::set(const Prepared& p, const Unit& u, const Option& o) noexcept {
        auto& ga = *p.pair.first.infection->latentGenotype(u.locus);
        auto& gb = *p.pair.second.infection->latentGenotype(u.locus);
        if (o.first == ga.value() and o.second == gb.value()) {
            return;
        }
        Proposal<T, Engine> move(*target_, *rng_, innerCounts_);
        if (o.first != ga.value()) {
            move.change(ga, o.first);
        }
        if (o.second != gb.value()) {
            move.change(gb, o.second);
        }
        move.accept();
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::pilot(const Prepared& p) noexcept {
        constexpr double least = 1e-9;
        std::vector<Option> options;
        std::vector<double> values;
        std::vector<double> gains(p.units.size());
        for (size_t u = 0; u < p.units.size(); ++u) {
            gains[u] = bestAt(p, p.units[u], options, values).first;
        }
        for (size_t step = 0; step < 4 * p.units.size() + 40; ++step) {
            const auto u = static_cast<size_t>(std::ranges::max_element(gains) - gains.begin());
            if (gains[u] <= least) {
                return;
            }
            const auto [gain, option] = bestAt(p, p.units[u], options, values);
            gains[u]                  = gain;
            if (gain <= least or std::ranges::any_of(gains, [&](const double other) { return other > gain; })) {
                continue;
            }
            set(p, p.units[u], option);
            for (size_t v = 0; v < p.units.size(); ++v) {
                if (p.units[v].locus == p.units[u].locus) {
                    gains[v] = bestAt(p, p.units[v], options, values).first;
                }
            }
        }
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    double GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::walk(const Prepared& p, const Genotypes* towards) noexcept {
        boost::random::uniform_01<> unit;
        std::vector<Option> options;
        std::vector<double> values;
        double logQ = 0;
        for (size_t u = 0; u < p.units.size(); ++u) {
            const auto& [locus, allele] = p.units[u];
            optionsAt(p, p.units[u], options, values);
            if (options.empty()) {
                return -std::numeric_limits<double>::infinity();
            }
            const double top = *std::ranges::max_element(values);
            if (not std::isfinite(top)) {
                return -std::numeric_limits<double>::infinity();
            }
            double total = 0;
            for (const double v : values) {
                total += std::exp(v - top);
            }
            size_t pick = 0;
            if (towards != nullptr) {
                const auto& [ta, tb] = (*towards)[u];
                pick                 = static_cast<size_t>(std::ranges::find_if(options, [&](const Option& o) { return o.first.allele(allele) == ta.allele(allele) and o.second.allele(allele) == tb.allele(allele); }) -
                                       options.begin());
                if (pick == options.size()) {
                    return -std::numeric_limits<double>::infinity();
                }
            } else {
                double r = unit(*rng_) * total;
                while (pick + 1 < options.size() and (r -= std::exp(values[pick] - top)) > 0) {
                    ++pick;
                }
            }
            logQ += values[pick] - top - std::log(total);
            set(p, p.units[u], options[pick]);
        }
        return logQ;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    typename GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::Genotypes GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::snapshot(const Prepared& p) const {
        Genotypes out;
        for (const auto& u : p.units) {
            out.emplace_back(p.pair.first.infection->latentGenotype(u.locus)->value(), p.pair.second.infection->latentGenotype(u.locus)->value());
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::restore(const Prepared& p, const double durationA, const double durationB, const Genotypes& genotypes) noexcept {
        const auto& a = *p.pair.first.infection;
        const auto& b = *p.pair.second.infection;
        Proposal<T, Engine> back(*target_, *rng_, innerCounts_);
        if (a.infectionDuration()->value() != durationA) {
            back.change(*a.infectionDuration(), durationA);
        }
        if (b.infectionDuration()->value() != durationB) {
            back.change(*b.infectionDuration(), durationB);
        }
        std::set<const void*> done;
        for (size_t u = 0; u < p.units.size(); ++u) {
            if (not done.insert(p.units[u].locus.get()).second) {
                continue;
            }
            if (a.latentGenotype(p.units[u].locus)->value() != genotypes[u].first) {
                back.change(*a.latentGenotype(p.units[u].locus), genotypes[u].first);
            }
            if (b.latentGenotype(p.units[u].locus)->value() != genotypes[u].second) {
                back.change(*b.latentGenotype(p.units[u].locus), genotypes[u].second);
            }
        }
        back.accept();
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::swap(const Prepared& p) noexcept {
        const auto& a = *p.pair.first.infection;
        const auto& b = *p.pair.second.infection;
        const double before    = target_->value();
        const double durationA = a.infectionDuration()->value(), durationB = b.infectionDuration()->value();
        const auto current     = snapshot(p);
        if (not exchange(p)) {
            countMove(false);
            return;
        }
        pilot(p);
        const double logForward = walk(p, nullptr);
        const double after      = target_->value();
        const double proposedA = a.infectionDuration()->value(), proposedB = b.infectionDuration()->value();
        const auto proposed    = snapshot(p);
        double logReverse      = -std::numeric_limits<double>::infinity();
        if (std::isfinite(logForward) and exchange(p)) {
            pilot(p);
            logReverse = walk(p, &current);
        }
        boost::random::uniform_01<> unit;
        const bool accepted = std::isfinite(logReverse) and std::log(unit(*rng_)) < after - before + logReverse - logForward;
        if (accepted) {
            restore(p, proposedA, proposedB, proposed);
        } else {
            restore(p, durationA, durationB, current);
        }
        countMove(accepted);
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_GUIDEDPAIRSAMPLER_H
