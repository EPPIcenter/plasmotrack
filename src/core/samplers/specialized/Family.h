//
// Families: the Mutual pairs whose Observed genotypes are close, and the connected triples of such pairs. The family
// moves (GuidedGroupSampler, PairLocusGibbsSampler) work on them. Built from the data only, so they are fixed.
//

#ifndef TRANSMISSION_NETWORKS_APP_FAMILY_H
#define TRANSMISSION_NETWORKS_APP_FAMILY_H

#include <algorithm>
#include <map>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    // Whether `child`'s Observed genotype at a Locus has an Allele `parent`'s lacks. Never at a Locus where either
    // has no Allele observed: a missing observation says nothing about the Alleles present.
    template<typename Genotype>
    bool novelAt(const Genotype& parent, const Genotype& child) {
        return parent.totalPositiveCount() > 0 and child.totalPositiveCount() > 0 and Genotype::falsePositiveCount(parent, child) > 0;
    }

    // The novelty of `child` over `parent`: the Loci, observed in both, where `child`'s Observed genotype has an Allele
    // `parent`'s lacks.
    template<typename InfectionEventImpl>
    unsigned int novelty(const InfectionEventImpl& parent, const InfectionEventImpl& child) {
        unsigned int n = 0;
        for (const auto& locus : child.loci()) {
            if (not parent.observedGenotype().contains(locus) or not child.observedGenotype().contains(locus)) continue;
            n += novelAt(parent.observedGenotype().at(locus)->value(), child.observedGenotype().at(locus)->value());
        }
        return n;
    }

    // Whether two Infections are close: the novelty in the direction where it is smaller is at most `maxNovelty`.
    template<typename InfectionEventImpl>
    bool close(const InfectionEventImpl& a, const InfectionEventImpl& b, const unsigned int maxNovelty) {
        return std::min(novelty(a, b), novelty(b, a)) <= maxNovelty;
    }

    // The close pairs among `pairs`, in their order. A Pair has members `first` and `second`, each with an `infection`.
    template<typename Pair>
    std::vector<Pair> closePairs(const std::vector<Pair>& pairs, const unsigned int maxNovelty) {
        std::vector<Pair> out;
        for (const auto& pair : pairs) {
            if (close(*pair.first.infection, *pair.second.infection, maxNovelty)) out.push_back(pair);
        }
        return out;
    }

    // The families of the close pairs `pairs`: each pair, then each connected triple (a member and two of its close
    // partners), in an order fixed by the members' uids.
    template<typename Pair>
    std::vector<std::vector<decltype(Pair::first)>> families(const std::vector<Pair>& pairs) {
        using Member = decltype(Pair::first);
        std::vector<std::vector<Member>> out;
        // Each close pair under both members, the member first.
        std::map<int, std::vector<std::pair<Member, Member>>> partnersOf;
        for (const auto& [a, b] : pairs) {
            out.push_back({a, b});
            partnersOf[a.infection->uid()].push_back({a, b});
            partnersOf[b.infection->uid()].push_back({b, a});
        }
        for (const auto& [centre, partners] : partnersOf) {
            for (size_t i = 0; i < partners.size(); ++i) {
                for (size_t j = i + 1; j < partners.size(); ++j) out.push_back({partners[i].first, partners[i].second, partners[j].second});
            }
        }
        return out;
    }

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_FAMILY_H
