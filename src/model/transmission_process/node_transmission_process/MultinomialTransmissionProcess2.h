//
// Created by Maxwell Murphy on 11/17/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"
#include "core/utils/ProbAnyMissing.h"
#include "core/utils/numerics.h"
#include <cmath>
#include <memory>

#include <boost/math/special_functions/binomial.hpp>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;
    using Probability = core::computation::Probability;

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    class MultinomialTransmissionProcess2 : public core::computation::Computation<std::array<Likelihood, MAX_STRAINS*(MAX_PARENTSET_SIZE + 1)>>,
                                            public core::abstract::Observable<MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>>,
                                            public core::abstract::Cacheable<MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>>,
                                            public core::abstract::Checkpointable<MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTSET_SIZE + 1)>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;

    public:
        // Transmission process that allows for mutations
        explicit MultinomialTransmissionProcess2(p_ParameterDouble mean_strains_transmitted);

        std::array<Likelihood, MAX_STRAINS*(MAX_PARENTSET_SIZE + 1)> value() override;

        Likelihood probNumStrains(int num_strains, int num_parents) {
            return this->value()[(num_parents - 1) * MAX_STRAINS + (num_strains - 1)];
        }

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet,
                                          std::shared_ptr<ParentSetSizePriorImpl> psp);


        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp,
                                          std::shared_ptr<ParentSetSizePriorImpl> psp);

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp,
                                          std::shared_ptr<ParentSetSizePriorImpl> psp);


    private:
        friend class core::abstract::Checkpointable<MultinomialTransmissionProcess2, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTSET_SIZE + 1)>>;
        friend class core::abstract::Cacheable<MultinomialTransmissionProcess2>;

        p_ParameterDouble mean_strains_transmitted_;
        core::utils::probAnyMissingFunctor probAnyMissing_;
        Probability mutationRate_ = 0.001;
    };

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::MultinomialTransmissionProcess2(p_ParameterDouble mean_strains_transmitted) : mean_strains_transmitted_(std::move(mean_strains_transmitted)) {
        mean_strains_transmitted_->registerCacheableCheckpointTarget(this);
        mean_strains_transmitted_->add_post_change_listener([this]() {
            this->setDirty();
        });

        this->value_.fill(-std::numeric_limits<Likelihood>::infinity());
        this->setDirty();
        this->MultinomialTransmissionProcess2::value();
    }

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizeLikelihood>
    std::array<Likelihood, MAX_STRAINS*(MAX_PARENTSET_SIZE + 1)> MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizeLikelihood>::value() {
        /*
         * Value is a matrix that gives the probability of transmitting some number of strains from some number of parents. The
         * rows are the number of parents, and the columns are the number of strains. We assume the number of strains transmitted
         * from a single parent is zt-poisson with rate equal to `mean_strains_transmitted_`. In the case of multiple parents, the number
         * of strains transmitted is the sum of the number of strains transmitted from each parent, description of the distribution may be found
         * in the paper https://www.semanticscholar.org/paper/On-the-sum-of-independent-zero-truncated-Poisson-Springael-Nieuwenhuyse/f3061cd20c1d1e765cf83b7978206c22bd94ce55#citing-papers
         * entitled "On the sum of independent zero-truncated Poisson random variables"
         * todo: mean_strains_transmitted isn't actually the mean of the distribution
         */
        if (this->isDirty()) {
            const double lambda = mean_strains_transmitted_->value();
            for (unsigned int kk = 0; kk < MAX_PARENTSET_SIZE + 1; ++kk) {
                const int num_parents = static_cast<int>(kk) + 1;
                Probability denominator = 0.0;
                for (unsigned int jj = kk; jj < MAX_STRAINS; ++jj) {
                    const int num_strains = static_cast<int>(jj) + 1;

                    double correction = 0;
                    for (int ii = 0; ii < num_parents; ++ii) {
                        correction += std::pow(-1, ii) * std::pow(num_parents - ii, num_strains) * boost::math::binomial_coefficient<double>(num_parents, ii);
                    }
                    this->value_[kk * MAX_STRAINS + jj] = num_strains * std::log(lambda) - (num_parents * std::log(std::exp(lambda) - 1)) - std::log(boost::math::factorial<double>(num_strains)) + std::log(correction);
                    denominator += exp(this->value_[kk * MAX_STRAINS + jj]);
                }

                // Normalize
                denominator = std::log(denominator);
                for (unsigned int jj = 0; jj < MAX_STRAINS; ++jj) {
                    this->value_[kk * MAX_STRAINS + jj] -= denominator;
                }
            }
            this->setClean();
        }
        return this->value_;
    }

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizeLikelihood>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizeLikelihood>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet, std::shared_ptr<ParentSetSizeLikelihood> psp) {
        const size_t numParents = parentSet.size();
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};


        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();
            std::vector<Probability> parent_pop_freqs(childGenotype.totalAlleles(), 0);

            for (const auto& parent : parentSet) {
                const auto& parentGenotype = parent->latentGenotype(locus)->value();
                // Every parent must have at least one allele in common with the child at each locus
                if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
                    return -std::numeric_limits<Likelihood>::infinity();
                }

                const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                    if (parentGenotype.allele(j)) {
                        parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / (totalAllelesPresent * numParents));
                    }
                }
            }

            long double total_positive = 0.0;
            for (size_t j = 0; j < childGenotype.totalAlleles(); ++j) {
                total_positive += parent_pop_freqs[j] > 0 ? 1.0 : 0.0;
            }
            const long double total_negative = childGenotype.totalAlleles() - total_positive;

            for (size_t j = 0; j < childGenotype.totalAlleles(); ++j) {
                if (parent_pop_freqs[j] == 0) {
                    parent_pop_freqs[j] = mutationRate_ / total_negative;
                } else {
                    parent_pop_freqs[j] -= mutationRate_ / total_positive;
                }
            }

            Probability constrainedSetProb = 0.0;
            std::vector<Likelihood> prVec{};
            prVec.reserve(childGenotype.totalAlleles());
            bool zeroProbEvent = false;

            for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
                if (childGenotype.allele(i) > 0) {
                    prVec.push_back(parent_pop_freqs[i]);
                    constrainedSetProb += parent_pop_freqs[i];
                    zeroProbEvent = zeroProbEvent || std::abs(parent_pop_freqs[i]) < 1e-10;
                }
            }

            if (zeroProbEvent) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            for (Likelihood& af : prVec) {
                af /= constrainedSetProb;
            }

            const Likelihood logConstrainedSetProb = std::log(constrainedSetProb);
            const std::vector<Likelihood>& pamVec = probAnyMissing_.vectorized(prVec, MAX_STRAINS);
            for (unsigned int numStrains = 1; numStrains <= MAX_STRAINS; ++numStrains) {
                unsigned int idx = numStrains - 1;

                if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                    continue;
                }
                logLikelihoods[idx] += pamVec[idx] >= 1.0 ? -std::numeric_limits<Likelihood>::infinity() : std::log(1.0 - pamVec[idx]) + logConstrainedSetProb * numStrains;
            }
        }

        Likelihood maxLl = -std::numeric_limits<Likelihood>::infinity();
        for (unsigned int numStrains = numParents; numStrains <= MAX_STRAINS; ++numStrains) {
            unsigned int idx = numStrains - 1;
            // Add the probability of the number of strains
            if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                continue;
            }
            logLikelihoods[idx] += this->probNumStrains(numStrains, numParents);
            maxLl = std::max(maxLl, logLikelihoods[idx]);
        }


        Likelihood llik = core::utils::logSumExpKnownMax(logLikelihoods.begin(), logLikelihoods.end(), maxLl);

        // Add the prior on the number of parents
        llik += psp->value()(numParents);
        return llik;
    }

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizeLikelihood>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizeLikelihood>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                                                                                                                                               std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                                                                                                                                               const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet,
                                                                                                                                                               std::shared_ptr<SourceTransmissionProcessImpl> stp,
                                                                                                                                                               std::shared_ptr<ParentSetSizeLikelihood> psp) {
        const size_t numParents = parentSet.size() + 1;// Add one for the latent parent
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};
        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();

            std::vector<Probability> parent_pop_freqs(childGenotype.totalAlleles(), 0.0);

            for (const auto& parent : parentSet) {
                const auto& parentGenotype = parent->latentGenotype(locus)->value();
                // Every parent must have at least one allele in common with the child at each locus
                if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
                    return -std::numeric_limits<Likelihood>::infinity();
                }

                const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                    // add an extra count to numParents to account for the latent parent
                    if (parentGenotype.allele(j)) {
                        parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / totalAllelesPresent / (numParents));
                    }
                }
            }


            const auto& latentParentGenotype = latentParent->latentGenotype(locus)->value();
            const int totalAllelesPresent = latentParentGenotype.totalPositiveCount();

            if (GeneticsImpl::truePositiveCount(latentParentGenotype, childGenotype) == 0) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            long double total_positive = 0.0;
            for (size_t j = 0; j < latentParentGenotype.totalAlleles(); ++j) {
                parent_pop_freqs[j] += (static_cast<Probability>(latentParentGenotype.allele(j)) / totalAllelesPresent / (numParents));
                total_positive += parent_pop_freqs[j] > 0 ? 1.0 : 0.0;
            }
            const long double total_negative = latentParentGenotype.totalAlleles() - total_positive;

            for (size_t j = 0; j < latentParentGenotype.totalAlleles(); ++j) {
                if (parent_pop_freqs[j] == 0) {
                    parent_pop_freqs[j] = mutationRate_ / total_negative;
                } else {
                    parent_pop_freqs[j] -= mutationRate_ / total_positive;
                }
            }

            Probability constrainedSetProb = 0.0;
            std::vector<Likelihood> prVec{};
            prVec.reserve(childGenotype.totalAlleles());
            bool zeroProbEvent = false;

            for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
                if (childGenotype.allele(i) > 0) {
                    prVec.push_back(parent_pop_freqs[i]);
                    constrainedSetProb += parent_pop_freqs[i];
                    zeroProbEvent = zeroProbEvent || std::abs(parent_pop_freqs[i]) < 1e-6;
                }
            }

            if (zeroProbEvent) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            for (Likelihood& af : prVec) {
                af /= constrainedSetProb;
            }

            const Likelihood logConstrainedSetProb = std::log(constrainedSetProb);
            const std::vector<Likelihood>& pamVec = probAnyMissing_.vectorized(prVec, MAX_STRAINS);
            for (unsigned int numStrains = 1; numStrains <= MAX_STRAINS; ++numStrains) {
                const unsigned int idx = numStrains - 1;

                if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                    continue;
                }
                logLikelihoods[idx] += pamVec[idx] >= 1.0 ? -std::numeric_limits<Likelihood>::infinity() : std::log(1.0 - pamVec[idx]) + logConstrainedSetProb * numStrains;
            }
        }

        Likelihood maxLl = -std::numeric_limits<Likelihood>::infinity();
        for (unsigned int numStrains = numParents; numStrains <= MAX_STRAINS; ++numStrains) {
            unsigned int idx = numStrains - 1;
            // Add the probability of the number of strains
            if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                continue;
            }
            logLikelihoods[idx] += this->probNumStrains(numStrains, numParents);
            maxLl = std::max(maxLl, logLikelihoods[idx]);
        }

        Likelihood llik = core::utils::logSumExpKnownMax(logLikelihoods.begin(), logLikelihoods.end(), maxLl) + stp->value();

        // Add the prior on the number of parents plus the source
        llik += psp->value()[numParents];
        return llik;
    }

    template<unsigned int MAX_PARENTSET_SIZE, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizeLikelihood>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess2<MAX_PARENTSET_SIZE, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizeLikelihood>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                                                                                                                                               std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                                                                                                                                               std::shared_ptr<SourceTransmissionProcessImpl> stp,
                                                                                                                                                               std::shared_ptr<ParentSetSizeLikelihood> psp) {
        constexpr size_t numParents = 1;
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};
        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();

            std::vector<Probability> parent_pop_freqs(childGenotype.totalAlleles(), 0.0);

            const auto& parentGenotype = latentParent->latentGenotype(locus)->value();
            // Every parent must have at least one allele in common with the child at each locus
            if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            const int totalAllelesPresent = parentGenotype.totalPositiveCount();
            long double total_positive = 0.0;
            for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                if (parentGenotype.allele(j)) {
                    parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / totalAllelesPresent);
                    total_positive += parent_pop_freqs[j] > 0 ? 1.0 : 0.0;
                }
            }

            const long double total_negative = parentGenotype.totalAlleles() - total_positive;

            for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                if (parent_pop_freqs[j] == 0) {
                    parent_pop_freqs[j] = mutationRate_ / total_negative;
                } else {
                    parent_pop_freqs[j] -= mutationRate_ / total_positive;
                }
            }

            Probability constrainedSetProb = 0.0;
            std::vector<Likelihood> prVec{};
            prVec.reserve(childGenotype.totalAlleles());
            bool zeroProbEvent = false;

            for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
                if (childGenotype.allele(i) > 0) {
                    prVec.push_back(parent_pop_freqs[i]);
                    constrainedSetProb += parent_pop_freqs[i];
                    zeroProbEvent = zeroProbEvent || std::abs(parent_pop_freqs[i]) < 1e-6;
                }
            }

            if (zeroProbEvent) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            for (Likelihood& af : prVec) {
                af /= constrainedSetProb;
            }

            const Likelihood logConstrainedSetProb = std::log(constrainedSetProb);
            const std::vector<Likelihood>& pamVec = probAnyMissing_.vectorized(prVec, MAX_STRAINS);
            for (unsigned int numStrains = 1; numStrains <= MAX_STRAINS; ++numStrains) {
                const unsigned int idx = numStrains - 1;
                if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                    continue;
                }
                logLikelihoods[idx] += pamVec[idx] >= 1.0 ? -std::numeric_limits<Likelihood>::infinity() : std::log(1.0 - pamVec[idx]) + logConstrainedSetProb * numStrains;
            }
        }

        Likelihood maxLl = -std::numeric_limits<Likelihood>::infinity();
        for (unsigned int numStrains = numParents; numStrains <= MAX_STRAINS; ++numStrains) {
            const unsigned int idx = numStrains - 1;
            // Add the probability of the number of strains
            if (logLikelihoods[idx] == -std::numeric_limits<Likelihood>::infinity()) {
                continue;
            }
            logLikelihoods[idx] += this->probNumStrains(numStrains, numParents);
            maxLl = std::max(maxLl, logLikelihoods[idx]);
        }

        Likelihood llik = core::utils::logSumExpKnownMax(logLikelihoods.begin(), logLikelihoods.end(), maxLl) + stp->value();

        // Add the prior on the number of parents plus the source
        llik += psp->value()[numParents];
        return llik;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H
