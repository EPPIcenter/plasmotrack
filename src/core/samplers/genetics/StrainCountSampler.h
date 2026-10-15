//
// Removes or adds one strain's alleles across an infection's loci in one step.
//

#ifndef TRANSMISSION_NETWORKS_APP_STRAINCOUNTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_STRAINCOUNTSAMPLER_H

#include "core/datatypes/Genotype.h"
#include "core/samplers/AbstractSampler.h"

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    /*
     * A Metropolis-Hastings move on an infection's latent genotypes at several loci at once, changing the largest
     * number of alleles at any locus (M, a lower bound on the infection's strains) by one.
     *
     * Remove (with probability 1/2, when M >= 2): at every locus with M alleles, one of them, chosen uniformly, is
     * removed. Afterwards the largest count is M - 1.
     *
     * Add (with probability 1/2): at each locus with M alleles, independently, one allele is added with probability
     * pObserved if the locus has an observed allele the genotype lacks and pOther otherwise. The allele is drawn from
     * the absent ones, an observed one weighted observedWeight times another. Afterwards the loci with M + 1 alleles are
     * exactly those an allele was added to. An add that chooses no locus changes nothing.
     *
     * Each is the other's reverse: a remove from the result of an add removes at exactly the loci the add chose, and an
     * add from the result of a remove can choose exactly the loci the remove changed. The Hastings ratio is the product
     * of those choices' probabilities. Genotypes with more than maxCOI alleles at a locus have probability zero and are
     * rejected before any state changes.
     *
     * An infection's strains are shared by its loci, so a strain the data don't need costs little at any one locus
     * once its alleles are present at the others. Removing it one locus at a time then earns almost nothing per step,
     * and the chain keeps it. This move removes, or adds, it at every locus at once.
     */
    template<typename T, typename Engine, typename InfectionEventImpl>
    class StrainCountSampler : public AbstractSampler {
        using Genotype = core::datatypes::Genotype;

    public:
        StrainCountSampler(std::shared_ptr<InfectionEventImpl> infection, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int maxCOI,
                           double pObserved = 0.5, double pOther = 0.01, double observedWeight = 10.0)
            : infection_(std::move(infection)), target_(std::move(target)), rng_(std::move(rng)), maxCOI_(maxCOI),
              pObserved_(pObserved), pOther_(pOther), observedWeight_(observedWeight) {}

        void update() noexcept override;

    private:
        static constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

        // The locus' observed alleles, or none.
        [[nodiscard]] Genotype observed(size_t l, const Genotype& like) const;
        // Log probability that an add from `genotypes` chooses exactly the loci in `chosen` (a flag per locus) and adds
        // allele `added[l]` at each chosen locus.
        [[nodiscard]] double logAddProb(const std::vector<Genotype>& genotypes, const std::vector<bool>& chosen, const std::vector<unsigned int>& added) const;
        // Log probability that a remove from `genotypes` removes allele `removed[l]` at each locus with the largest count.
        [[nodiscard]] static double logRemoveProb(const std::vector<Genotype>& genotypes);
        [[nodiscard]] static unsigned int largest(const std::vector<Genotype>& genotypes);

        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int maxCOI_;
        double pObserved_;
        double pOther_;
        double observedWeight_;

        boost::random::uniform_01<> uniform_{};
    };

    template<typename T, typename Engine, typename InfectionEventImpl>
    typename StrainCountSampler<T, Engine, InfectionEventImpl>::Genotype StrainCountSampler<T, Engine, InfectionEventImpl>::observed(const size_t l, const Genotype& like) const {
        const auto& locus = infection_->loci()[l];
        if (infection_->observedGenotype().contains(locus)) {
            return infection_->observedGenotype(locus)->value();
        }
        auto none = like;
        none.reset();
        return none;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    unsigned int StrainCountSampler<T, Engine, InfectionEventImpl>::largest(const std::vector<Genotype>& genotypes) {
        unsigned int m = 0;
        for (const auto& g : genotypes) {
            m = std::max(m, g.totalPositiveCount());
        }
        return m;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double StrainCountSampler<T, Engine, InfectionEventImpl>::logRemoveProb(const std::vector<Genotype>& genotypes) {
        const unsigned int m = largest(genotypes);
        double out           = 0;
        for (const auto& g : genotypes) {
            if (g.totalPositiveCount() == m) {
                out -= std::log(static_cast<double>(m));
            }
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    double StrainCountSampler<T, Engine, InfectionEventImpl>::logAddProb(const std::vector<Genotype>& genotypes, const std::vector<bool>& chosen, const std::vector<unsigned int>& added) const {
        const unsigned int m = largest(genotypes);
        double out           = 0;
        for (size_t l = 0; l < genotypes.size(); ++l) {
            const auto& g = genotypes[l];
            if (g.totalPositiveCount() != m) {
                if (chosen[l]) {
                    return NEG_INF;
                }
                continue;
            }
            const auto obs  = observed(l, g);
            double observedAbsent = 0, otherAbsent = 0;
            for (unsigned int a = 0; a < g.totalAlleles(); ++a) {
                if (not g.allele(a)) {
                    (obs.allele(a) ? observedAbsent : otherAbsent) += 1;
                }
            }
            const double p = observedAbsent > 0 ? pObserved_ : pOther_;
            if (not chosen[l]) {
                out += std::log1p(-p);
                continue;
            }
            const double total = observedWeight_ * observedAbsent + otherAbsent;
            if (total <= 0 or g.allele(added[l])) {
                return NEG_INF;
            }
            const double weight = obs.allele(added[l]) ? observedWeight_ : 1.0;
            out += std::log(p) + std::log(weight / total);
        }
        return out;
    }

    template<typename T, typename Engine, typename InfectionEventImpl>
    void StrainCountSampler<T, Engine, InfectionEventImpl>::update() noexcept {
        const auto& loci = infection_->loci();
        std::vector<Genotype> current;
        current.reserve(loci.size());
        for (const auto& locus : loci) {
            current.push_back(infection_->latentGenotype(locus)->value());
        }
        const unsigned int m = largest(current);
        auto proposed        = current;
        std::vector<bool> chosen(loci.size(), false);
        std::vector<unsigned int> changed(loci.size(), 0);// the allele removed or added at each chosen locus
        double logForward = 0, logReverse = 0;

        if (uniform_(*rng_) < 0.5) {
            // Remove one allele at every locus with m alleles.
            if (m < 2) {
                return;
            }
            for (size_t l = 0; l < loci.size(); ++l) {
                if (current[l].totalPositiveCount() != m) {
                    continue;
                }
                std::vector<unsigned int> present;
                for (unsigned int a = 0; a < current[l].totalAlleles(); ++a) {
                    if (current[l].allele(a)) {
                        present.push_back(a);
                    }
                }
                boost::random::uniform_int_distribution<size_t> pick(0, present.size() - 1);
                changed[l] = present[pick(*rng_)];
                chosen[l]  = true;
                proposed[l].set(changed[l], false);
            }
            logForward = logRemoveProb(current);
            logReverse = logAddProb(proposed, chosen, changed);
        } else {
            // Add one allele at each chosen locus with m alleles.
            if (m + 1 > maxCOI_) {
                return;
            }
            bool any = false;
            for (size_t l = 0; l < loci.size(); ++l) {
                const auto& g = current[l];
                if (g.totalPositiveCount() != m) {
                    continue;
                }
                const auto obs = observed(l, g);
                std::vector<unsigned int> absent;
                std::vector<double> weights;
                for (unsigned int a = 0; a < g.totalAlleles(); ++a) {
                    if (not g.allele(a)) {
                        absent.push_back(a);
                        weights.push_back(obs.allele(a) ? observedWeight_ : 1.0);
                    }
                }
                if (absent.empty()) {
                    continue;
                }
                const bool observedAbsent = std::ranges::any_of(absent, [&](const unsigned int a) { return obs.allele(a); });
                if (uniform_(*rng_) >= (observedAbsent ? pObserved_ : pOther_)) {
                    continue;
                }
                double total = 0;
                for (const double w : weights) {
                    total += w;
                }
                double u   = uniform_(*rng_) * total;
                size_t pick = absent.size() - 1;
                for (size_t i = 0; i < absent.size(); ++i) {
                    u -= weights[i];
                    if (u < 0) {
                        pick = i;
                        break;
                    }
                }
                changed[l] = absent[pick];
                chosen[l]  = true;
                proposed[l].set(changed[l], true);
                any = true;
            }
            if (not any) {
                return;// nothing changes
            }
            logForward = logAddProb(current, chosen, changed);
            logReverse = logRemoveProb(proposed);
        }

        auto proposal = this->propose(*target_, *rng_);
        for (size_t l = 0; l < loci.size(); ++l) {
            if (chosen[l]) {
                proposal.change(*infection_->latentGenotype(loci[l]), proposed[l]);
            }
        }
        proposal.metropolis(logReverse - logForward);
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_STRAINCOUNTSAMPLER_H
