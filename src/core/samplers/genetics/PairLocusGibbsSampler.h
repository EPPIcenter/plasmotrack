//
// Changes two related infections' latent genotypes at one locus together, by a Gibbs step over a fixed
// set of candidates built from the observed data.
//

#ifndef TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/Proposal.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    /*
     * Some traps need two infections' genotypes at one locus to change at once: one infection's observed allele
     * called a false positive while it takes its relative's allele, and the relative dropping a hidden allele it held
     * for the first. Each change alone costs tens of log units.
     *
     * Each (pair, locus) unit has, per member, a fixed candidate set from the observed data only: the member's
     * observed genotype, the partner's observed genotype, their union and their intersection (when not empty), those
     * within the COI bound; a member unobserved at the locus also every single-allele genotype from the partner's
     * observed alleles. An update picks a unit uniformly; if the members' current genotypes are both in their
     * candidate sets, it draws the pair's genotypes from the target restricted to the product of the sets (a Gibbs step
     * on a fixed set containing the current state, which leaves the target invariant); otherwise it does nothing.
     */
    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    class PairLocusGibbsSampler : public AbstractSampler {
    public:
        using Infection = std::shared_ptr<InfectionEventImpl>;

        PairLocusGibbsSampler(const std::vector<std::pair<Infection, Infection>>& pairs, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int maxCOI, int moves = 1, unsigned int maxNovelty = 3)
            : target_(std::move(target)), rng_(std::move(rng)), moves_(moves) {
            for (const auto& [a, b] : pairs) {
                unsigned int nab = 0, nba = 0;
                for (const auto& [locus, oa] : a->observedGenotype()) {
                    if (b->observedGenotype().contains(locus)) {
                        const auto& ob = b->observedGenotype().at(locus)->value();
                        nab += Genotype::falsePositiveCount(oa->value(), ob) > 0;
                        nba += Genotype::falsePositiveCount(ob, oa->value()) > 0;
                    }
                }
                if (std::min(nab, nba) > maxNovelty) continue;

                for (const auto& locus : a->loci()) {
                    const bool ha = a->observedGenotype().contains(locus) and a->observedGenotype().at(locus)->value().totalPositiveCount() > 0;
                    const bool hb = b->observedGenotype().contains(locus) and b->observedGenotype().at(locus)->value().totalPositiveCount() > 0;
                    if (not ha and not hb) continue;
                    if (ha and hb) {
                        const auto& oa = a->observedGenotype().at(locus)->value();
                        const auto& ob = b->observedGenotype().at(locus)->value();
                        if (oa == ob) continue;// nothing to exchange
                        const bool conflictAB = (nab <= maxNovelty) and (Genotype::falsePositiveCount(oa, ob) > 0);
                        const bool conflictBA = (nba <= maxNovelty) and (Genotype::falsePositiveCount(ob, oa) > 0);
                        if (not conflictAB and not conflictBA) continue;
                    }
                    Unit u{.a = a, .b = b, .locus = locus};
                    std::vector<Genotype> base;
                    if (ha) base.push_back(a->observedGenotype().at(locus)->value());
                    if (hb) base.push_back(b->observedGenotype().at(locus)->value());
                    if (ha and hb) {
                        base.push_back(Genotype::any(base[0], base[1]));
                        if (const auto s = Genotype::shared(base[0], base[1]); s.totalPositiveCount() > 0) base.push_back(s);
                    }
                    for (const auto& g : base) {
                        if (g.totalPositiveCount() == 0 or g.totalPositiveCount() > maxCOI) continue;
                        if (std::ranges::find(u.ca, g) == u.ca.end()) u.ca.push_back(g);
                        if (std::ranges::find(u.cb, g) == u.cb.end()) u.cb.push_back(g);
                    }
                    // A member unobserved at the locus may also hold any single allele observed in its partner there.
                    if (not ha or not hb) {
                        const auto& any = ha ? a->observedGenotype().at(locus)->value() : b->observedGenotype().at(locus)->value();
                        for (unsigned int allele = 0; allele < any.totalAlleles(); ++allele) {
                            if (not any.allele(allele)) continue;
                            Genotype single = any;
                            single.reset();
                            single.set(allele);
                            if (not ha and std::ranges::find(u.ca, single) == u.ca.end()) u.ca.push_back(single);
                            if (not hb and std::ranges::find(u.cb, single) == u.cb.end()) u.cb.push_back(single);
                        }
                    }
                    if (u.ca.size() * u.cb.size() > 1) {
                        units_.push_back(std::move(u));
                    }
                }
            }
        }

        void update() noexcept override {
            for (int m = 0; m < moves_ and not units_.empty(); ++m) {
                const auto& u = units_[boost::random::uniform_int_distribution<size_t>(0, units_.size() - 1)(*rng_)];
                auto& pa      = *u.a->latentGenotype(u.locus);
                auto& pb      = *u.b->latentGenotype(u.locus);
                const auto ca = pa.value(), cb = pb.value();
                if (std::ranges::find(u.ca, ca) == u.ca.end() or std::ranges::find(u.cb, cb) == u.cb.end()) {
                    countProposal(false);
                    continue;
                }
                std::vector<double> values;
                values.reserve(u.ca.size() * u.cb.size());
                for (const auto& ga : u.ca) {
                    for (const auto& gb : u.cb) {
                        Proposal<T, Engine> move(*target_, *rng_, inner_);
                        if (ga != ca) move.change(pa, ga);
                        if (gb != cb) move.change(pb, gb);
                        values.push_back(move.value());
                        move.revert();
                    }
                }
                const double top = *std::ranges::max_element(values);
                double total     = 0;
                for (const double v : values) total += std::exp(v - top);
                boost::random::uniform_01<> unit;
                double r    = unit(*rng_) * total;
                size_t pick = 0;
                while (pick + 1 < values.size() and (r -= std::exp(values[pick] - top)) > 0) ++pick;
                const auto& ga = u.ca[pick / u.cb.size()];
                const auto& gb = u.cb[pick % u.cb.size()];
                if (ga == ca and gb == cb) {
                    countProposal(false);
                    continue;
                }
                Proposal<T, Engine> move(*target_, *rng_, inner_);
                if (ga != ca) move.change(pa, ga);
                if (gb != cb) move.change(pb, gb);
                move.value();
                move.accept();
                countProposal(true);
            }
        }

        [[nodiscard]] size_t units() const noexcept { return units_.size(); }

    private:
        using Genotype = std::remove_cvref_t<decltype(std::declval<InfectionEventImpl&>().latentGenotype(std::declval<std::shared_ptr<LocusImpl>>())->value())>;
        struct Unit {
            Infection a, b;
            std::shared_ptr<LocusImpl> locus;
            std::vector<Genotype> ca, cb;
        };
        std::vector<Unit> units_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        int moves_;
        SamplerCounts inner_{};
    };

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H
