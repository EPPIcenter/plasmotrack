//
// Changes a close Mutual pair's Latent genotypes at one Locus together, by a Gibbs Proposal over a fixed set of
// candidates built from the Observed genotypes.
//

#ifndef TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/LogWeightedDraw.h"
#include "core/samplers/Proposal.h"
#include "core/samplers/specialized/Family.h"
#include "core/samplers/specialized/GuidedPairSampler.h"

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
     * called a false positive while it takes its relative's allele, and the relative dropping an allele it held unobserved
     * for the first. Each change alone costs tens of log units.
     *
     * Each (pair, locus) unit has, per member, a fixed candidate set from the observed data only: the member's
     * observed genotype, the partner's observed genotype, their union and their intersection (when not empty), those
     * within the COI bound; a member unobserved at the locus also every single-allele genotype from the partner's
     * observed alleles. An update picks a unit uniformly; if the members' current genotypes are both in their
     * candidate sets, it draws the pair's genotypes from the target restricted to the product of the sets (a Gibbs
     * Proposal on a fixed set containing the current state, which leaves the target invariant); otherwise it does
     * nothing.
     *
     * Only close pairs (specialized::close()) get units, and at a Locus both observe differently only when the
     * Observed genotypes there count toward a direction's novelty that is at most `maxNovelty`.
     */
    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    class PairLocusGibbsSampler : public AbstractSampler {
    public:
        using Infection = std::shared_ptr<InfectionEventImpl>;
        using Pair      = typename specialized::GuidedPairSampler<T, Engine, InfectionEventImpl, LocusImpl>::Pair;

        PairLocusGibbsSampler(const std::vector<Pair>& pairs, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int maxCOI, int proposals, unsigned int maxNovelty)
            : target_(std::move(target)), rng_(std::move(rng)), proposals_(proposals) {
            for (const auto& pair : pairs) {
                const auto& a = pair.first.infection;
                const auto& b = pair.second.infection;
                if (not specialized::close(*a, *b, maxNovelty)) continue;
                const unsigned int nab = specialized::novelty(*a, *b), nba = specialized::novelty(*b, *a);

                for (const auto& locus : a->loci()) {
                    const bool ha = a->observedGenotype().contains(locus) and a->observedGenotype().at(locus)->value().totalPositiveCount() > 0;
                    const bool hb = b->observedGenotype().contains(locus) and b->observedGenotype().at(locus)->value().totalPositiveCount() > 0;
                    if (not ha and not hb) continue;
                    if (ha and hb) {
                        const auto& oa = a->observedGenotype().at(locus)->value();
                        const auto& ob = b->observedGenotype().at(locus)->value();
                        if (oa == ob) continue;// nothing to exchange
                        const bool conflictAB = (nab <= maxNovelty) and specialized::novelAt(oa, ob);
                        const bool conflictBA = (nba <= maxNovelty) and specialized::novelAt(ob, oa);
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
                        if (std::ranges::find(u.candidatesA, g) == u.candidatesA.end()) u.candidatesA.push_back(g);
                        if (std::ranges::find(u.candidatesB, g) == u.candidatesB.end()) u.candidatesB.push_back(g);
                    }
                    // A member unobserved at the locus may also hold any single allele observed in its partner there.
                    if (not ha or not hb) {
                        const auto& any = ha ? a->observedGenotype().at(locus)->value() : b->observedGenotype().at(locus)->value();
                        for (unsigned int allele = 0; allele < any.totalAlleles(); ++allele) {
                            if (not any.allele(allele)) continue;
                            Genotype single = any;
                            single.reset();
                            single.set(allele);
                            if (not ha and std::ranges::find(u.candidatesA, single) == u.candidatesA.end()) u.candidatesA.push_back(single);
                            if (not hb and std::ranges::find(u.candidatesB, single) == u.candidatesB.end()) u.candidatesB.push_back(single);
                        }
                    }
                    if (u.candidatesA.size() * u.candidatesB.size() > 1) {
                        units_.push_back(std::move(u));
                    }
                }
            }
        }

        void update() noexcept override {
            for (int m = 0; m < proposals_ and not units_.empty(); ++m) {
                const auto& u = units_[boost::random::uniform_int_distribution<size_t>(0, units_.size() - 1)(*rng_)];
                auto& pa      = *u.a->latentGenotype(u.locus);
                auto& pb      = *u.b->latentGenotype(u.locus);
                const auto currentA = pa.value(), currentB = pb.value();
                if (std::ranges::find(u.candidatesA, currentA) == u.candidatesA.end() or std::ranges::find(u.candidatesB, currentB) == u.candidatesB.end()) {
                    countProposal(false);
                    continue;
                }
                const double currentVal = target_->value();
                std::vector<double> values;
                values.reserve(u.candidatesA.size() * u.candidatesB.size());
                for (const auto& ga : u.candidatesA) {
                    for (const auto& gb : u.candidatesB) {
                        if (ga == currentA and gb == currentB) {
                            values.push_back(currentVal);
                        } else {
                            Proposal<T, Engine> move(*target_, *rng_, inner_);
                            if (ga != currentA) move.change(pa, ga);
                            if (gb != currentB) move.change(pb, gb);
                            values.push_back(move.value());
                            move.revert();
                        }
                    }
                }
                const size_t pick = drawLogWeighted(values, logWeights(values), *rng_);
                const auto& ga = u.candidatesA[pick / u.candidatesB.size()];
                const auto& gb = u.candidatesB[pick % u.candidatesB.size()];
                if (ga == currentA and gb == currentB) {
                    countProposal(false);
                    continue;
                }
                Proposal<T, Engine> move(*target_, *rng_, inner_);
                if (ga != currentA) move.change(pa, ga);
                if (gb != currentB) move.change(pb, gb);
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
            std::vector<Genotype> candidatesA, candidatesB;
        };
        std::vector<Unit> units_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        int proposals_;
        SamplerCounts inner_{};
    };

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_PAIRLOCUSGIBBSSAMPLER_H
