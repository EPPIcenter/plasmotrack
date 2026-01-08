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

#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <memory>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;
    using Probability = core::computation::Probability;

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    class MultinomialTransmissionProcess : public core::computation::Computation<std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>>,
                                           public core::abstract::Observable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>>,
                                           public core::abstract::Cacheable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>>,
                                           public core::abstract::Checkpointable<MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
        using p_ParameterArray = std::shared_ptr<core::parameters::Parameter<std::array<Probability, MAX_PARENTS + 1>>>;

        template<typename GeneticsImpl>
        using p_Infection = std::shared_ptr<core::containers::Infection<GeneticsImpl>>;

        template<typename GeneticsImpl>
        using ParentSet = core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>;

        using p_SourceTransmissionProcess = std::shared_ptr<SourceTransmissionProcessImpl>;
        using p_ParentSetSizePrior = std::shared_ptr<ParentSetSizePriorImpl>;

    public:
        explicit MultinomialTransmissionProcess(p_ParameterDouble mean_strains_transmitted);

        std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)> value() override;

        Likelihood probNumStrains(const int num_strains, const int num_parents) {
            return this->value()[(num_parents - 1) * MAX_STRAINS + (num_strains - 1)];
        }

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(p_Infection<GeneticsImpl> infection, const ParentSet<GeneticsImpl>& parentSet, p_ParentSetSizePrior psp);


        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(p_Infection<GeneticsImpl> infection, p_Infection<GeneticsImpl> latentParent, const ParentSet<GeneticsImpl>& parentSet, p_SourceTransmissionProcess stp, p_ParentSetSizePrior psp);

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(p_Infection<GeneticsImpl> infection, p_Infection<GeneticsImpl> latentParent, p_SourceTransmissionProcess stp, p_ParentSetSizePrior psp);


    private:
        friend class core::abstract::Checkpointable<MultinomialTransmissionProcess, std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)>>;
        friend class core::abstract::Cacheable<MultinomialTransmissionProcess>;

        p_ParameterDouble meanStrainsTransmitted_;
        p_ParameterArray probParentSetSize_;

        core::utils::probAnyMissingFunctor probAnyMissing_;
    };

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::MultinomialTransmissionProcess(p_ParameterDouble mean_strains_transmitted) : meanStrainsTransmitted_(std::move(mean_strains_transmitted)) {
        meanStrainsTransmitted_->registerCacheableCheckpointTarget(this);
        meanStrainsTransmitted_->add_post_change_listener([this]() {
            this->setDirty();
        });

        this->value_.fill(-std::numeric_limits<Likelihood>::infinity());
        this->setDirty();
        this->MultinomialTransmissionProcess::value();
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    std::array<Likelihood, MAX_STRAINS*(MAX_PARENTS + 1)> MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::value() {
        /*
         * Value is a matrix that gives the probability of transmitting some number of strains from some number of parents. The
         * rows are the number of parents, and the columns are the number of strains. We assume the number of strains transmitted
         * from a single parent is zt-poisson with rate equal to `meanStrainsTransmitted_`. In the case of multiple parents, the number
         * of strains transmitted is the sum of the number of strains transmitted from each parent
         * todo: meanStrainsTransmitted_ isn't actually the mean of the distribution
         */
        if (this->isDirty()) {
            const double lambda = meanStrainsTransmitted_->value();
            for (unsigned int kk = 0; kk < MAX_PARENTS + 1; ++kk) {
                const int num_parents = static_cast<int>(kk) + 1;
                Probability denominator = 0.0;
                for (unsigned int jj = kk; jj < MAX_STRAINS; ++jj) {
                    const int num_strains = static_cast<int>(jj) + 1;

                    double correction = 0;
                    for (int ii = 0; ii < num_parents; ++ii) {
                        correction += std::pow(-1, ii) * std::pow(num_parents - ii, num_strains) * boost::math::binomial_coefficient<double>(num_parents, ii);
                    }
                    this->value_[kk * MAX_STRAINS + jj] = num_strains * std::log(lambda) - (num_parents * std::log(std::exp(lambda) - 1)) - std::log(boost::math::factorial<double>(num_strains)) + std::log(correction);
                    denominator += std::exp(this->value_[kk * MAX_STRAINS + jj]);
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



    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::calculateLogLikelihood(p_Infection<GeneticsImpl> infection, const ParentSet<GeneticsImpl>& parentSet, p_ParentSetSizePrior psp) {
        const size_t numParents = parentSet.size();
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};


        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();
            std::vector<Probability> parent_pop_freqs(childGenotype.totalAlleles(), 0);


            for (const auto& parent : parentSet) {
                const auto& parentGenotype = parent->latentGenotype(locus)->value();
                // Every parent must have at least one allele in common with the child at each locus
                // if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
                //     return -std::numeric_limits<Likelihood>::infinity();
                // }

                const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                    if (parentGenotype.allele(j)) {
                        parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / (totalAllelesPresent * numParents));
                    }
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

            if (prVec.empty() || zeroProbEvent) {
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

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::calculateLogLikelihood(
        p_Infection<GeneticsImpl> infection,
        p_Infection<GeneticsImpl> latentParent,
        const ParentSet<GeneticsImpl>& parentSet,
        p_SourceTransmissionProcess stp,
        p_ParentSetSizePrior psp) {
        const size_t numParents = parentSet.size() + 1;// Add one for the latent parent
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};
        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();

            std::vector<Likelihood> parent_pop_freqs(childGenotype.totalAlleles(), 0.0);

            for (const auto& parent : parentSet) {
                const auto& parentGenotype = parent->latentGenotype(locus)->value();
                // Every parent must have at least one allele in common with the child at each locus
                // if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
                //     return -std::numeric_limits<Likelihood>::infinity();
                // }


                const int totalAllelesPresent = parentGenotype.totalPositiveCount();
                for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                    if (parentGenotype.allele(j)) {
                        parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / totalAllelesPresent / static_cast<Probability>(numParents));
                    }
                }
            }


            const auto& latentParentGenotype = latentParent->latentGenotype(locus)->value();
            const int totalAllelesPresent = latentParentGenotype.totalPositiveCount();

            if (GeneticsImpl::truePositiveCount(latentParentGenotype, childGenotype) == 0) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            for (size_t j = 0; j < latentParentGenotype.totalAlleles(); ++j) {
                parent_pop_freqs[j] += (static_cast<Probability>(latentParentGenotype.allele(j)) / totalAllelesPresent / static_cast<Probability>(numParents));
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

            if (prVec.empty() || zeroProbEvent) {
                // print an error message
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

        // Add the prior on the number of parents
        llik += psp->value()(numParents);
        return llik;
    }

    template<unsigned int MAX_PARENTS, unsigned int MAX_STRAINS, typename SourceTransmissionProcessImpl, typename ParentSetSizePriorImpl>
    template<typename GeneticsImpl>
    Likelihood MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionProcessImpl, ParentSetSizePriorImpl>::calculateLogLikelihood(
            p_Infection<GeneticsImpl> infection,
            p_Infection<GeneticsImpl> latentParent,
            p_SourceTransmissionProcess stp,
            p_ParentSetSizePrior psp
            ) {

        constexpr size_t numParents = 1;
        const auto& loci = infection->loci();

        std::array<Likelihood, MAX_STRAINS> logLikelihoods{0};
        for (const auto& locus : loci) {
            const auto& childGenotype = infection->latentGenotype(locus)->value();

            std::vector<Probability> parent_pop_freqs(childGenotype.totalAlleles(), 0.0);

            const auto& parent = latentParent;
            const auto& parentGenotype = parent->latentGenotype(locus)->value();
            // Every parent must have at least one allele in common with the child at each locus
            // if (GeneticsImpl::truePositiveCount(parentGenotype, childGenotype) == 0) {
            //     return -std::numeric_limits<Likelihood>::infinity();
            // }

            const int totalAllelesPresent = parentGenotype.totalPositiveCount();
            for (size_t j = 0; j < parentGenotype.totalAlleles(); ++j) {
                if (parentGenotype.allele(j)) {
                    parent_pop_freqs[j] += (static_cast<Probability>(parentGenotype.allele(j)) / totalAllelesPresent);
                }
            }

            Probability constrainedSetProb = 0.0;
            std::vector<Likelihood> prVec{};
            prVec.reserve(childGenotype.totalAlleles());
            bool zeroProbEvent = false;

            for (size_t i = 0; i < childGenotype.totalAlleles(); ++i) {
                if (childGenotype.allele(i)) {
                    prVec.push_back(parent_pop_freqs[i]);
                    constrainedSetProb += parent_pop_freqs[i];
                    zeroProbEvent = zeroProbEvent || std::abs(parent_pop_freqs[i]) < 1e-6;
                }
            }

            if (prVec.empty() || zeroProbEvent) {
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

        // Add the prior on the number of parents
        llik += psp->value()(numParents);
        return llik;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_MULTINOMIALTRANSMISSIONPROCESS_H
