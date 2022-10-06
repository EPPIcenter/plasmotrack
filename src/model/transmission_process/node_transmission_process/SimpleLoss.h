//
// Created by Maxwell Murphy on 7/17/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H
#define TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"

#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"

#include "model/transmission_process/NetworkBasedTransmissionProcess.h"

#include <memory>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    class SimpleLoss : public core::computation::Computation<std::array<long double, MAX_TRANSMISSIONS + 1>>,
                       public core::abstract::Observable<SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>,
                       public core::abstract::Cacheable<SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>,
                       public core::abstract::Checkpointable<SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>, std::array<long double, MAX_TRANSMISSIONS + 1>> {

        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;

    public:
        explicit SimpleLoss(p_ParameterDouble loss_prob, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb);
//        SimpleLoss();

        std::array<long double, MAX_TRANSMISSIONS + 1> value() noexcept override;

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet);


        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp);

        template<typename GeneticsImpl>
        Likelihood calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp);


        template<typename GeneticsImpl>
        Likelihood peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                              const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet);


        template<typename GeneticsImpl>
        Likelihood peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp);

        template<typename GeneticsImpl>
        Likelihood peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection,
                                          std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent,
                                          std::shared_ptr<SourceTransmissionProcessImpl> stp);


    private:
        friend class core::abstract::Checkpointable<SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>, std::array<long double, MAX_TRANSMISSIONS + 1>>;
        friend class core::abstract::Cacheable<SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>;


        std::array<std::array<unsigned int, MAX_PARENTSET_SIZE + 1>, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> kVecs_ = core::utils::initKvecs<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE>();

        p_ParameterDouble lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb_;

        // Private store for probabilities when calculating the likelihood.
        std::array<long double, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> probs_{};

        // tracks [x, y, z] where x, y, and z are the number of alleles that are lost in the transmission process.
        std::array<unsigned int, MAX_PARENTSET_SIZE + 1> allelesLostCounter_{0};

        // For example, parentset size = 3, tracks [a, b, c, ab, bc, ac, abc] where a, b, c, ab, bc, ac, and abc are the number of alleles that are lost in the transmission process.
        std::array<unsigned int, core::utils::const_pow(2, MAX_PARENTSET_SIZE + 1) - 1> jointAllelesLostCounter_{0};
    };


    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::SimpleLoss(p_ParameterDouble loss_prob, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb) : lossProb_(std::move(loss_prob)), interTransmissionProb_(std::move(interTransmissionProb)) {
        lossProb_->template registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        interTransmissionProb_->registerCacheableCheckpointTarget(this);
        interTransmissionProb_->add_set_dirty_listener([this]() {
            this->setDirty();
        });

        this->value_.fill(0.0);
        this->setDirty();
        this->value();
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    std::array<long double, MAX_TRANSMISSIONS + 1> SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::value() noexcept {
        /*
         * Value is probability vector of not losing an allele up to that number of transmissions. Probability of losing an allele is 1 - value().
         */
        if (this->isDirty()) {
            this->value_[0] = 1.0;
            for (unsigned int i = 1; i <= MAX_TRANSMISSIONS; i++) {
                this->value_[i] = this->value_[i - 1] * (1.0 - this->lossProb_->value());
            }
            this->setClean();
        }
        return this->value_;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE> parentGenotypes{};


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;

            auto mutationFlag = childGenotype;

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                mutationFlag = mutationFlag.mutationMask(parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            if (mutationFlag.totalPositiveCount() > 0) {
//                fmt::print("Mutation flag is not zero: {}\n", core::io::serialize(mutationFlag));
                return -std::numeric_limits<Likelihood>::infinity();
            }


            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents; ++parentIdx) {
                psIdxGen.reset(numParents, parentIdx);
                while (!psIdxGen.completed) {
                    base = childGenotype;

                    // Invert parent genotypes by the index of the selected subset
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // AND all (maybe inverted) parent genotypes with the inverted child genotype
                    for (size_t ii = 0; ii < numParents; ++ii) {
                        base = GeneticsImpl::shared(base, parentGenotypes[ii]);
                    }

                    // Record the number of alleles that would be lost under that transmission
                    jointAllelesLostCounter_[jointIdx] += base.totalPositiveCount();
                    ++jointIdx;

                    // Invert again to get back to original state
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // Go to the next combination of parent indices
                    psIdxGen.next();
                }
            }
        }

        for (size_t ii = 0; ii < (unsigned long) std::pow(MAX_TRANSMISSIONS, numParents); ++ii) {
            auto kVec = kVecs_[ii];

            double all_parents_lost_prob = 1.0;
            // Calculate the probability of losses in the transmission process
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = 1.0 - this->value()[kVec[parentIdx]];
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->value()(kVec[parentIdx]));
            }

            // Calculate the probability of retentions in the transmission process
            probs_[ii] += jointAllelesLostCounter_[0] * std::log(1 - all_parents_lost_prob);
            size_t jointEventCounter = 1;
            for (size_t parentCount = 1; parentCount < numParents; ++parentCount) {

                // Now we need to iterate over the subsets of events
                psIdxGen.reset(numParents, parentCount);
                while (!psIdxGen.completed) {
                    // start with everything being multiplied out, i.e. all events happened
                    // for example (abc)
                    double tmpProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset
                        tmpProb /= (1.0 - this->value()[kVec[pIdx]]);
                    }
                    // what were left is the probability of the subset of parents retaining alleles
                    // i.e. (ab) if c was removed
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);

//        fmt::print("llik: {}\n", llik);
        if (std::isnan(llik) or llik <= -std::numeric_limits<Likelihood>::infinity()) {
//            fmt::print("NaN Encountered in SimpleLoss::calculateLogLikelihood without parent: {}\n", llik);
//            fmt::print("allelesLost: {}\n", core::io::serialize(allelesLostCounter_));
//            fmt::print("loss: {}\n", this->value()[1]);
//            fmt::print("retained: {}\n", 1.0 - this->value()[1]);
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

//        if (std::isnan(llik)) {
//            llik = -std::numeric_limits<Likelihood>::infinity();
//        }

        return llik;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE + 1> parentGenotypes{};



        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;

            auto mutationFlag = childGenotype;

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                mutationFlag = mutationFlag.mutationMask(parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            // calculate the latentParent contrib
            parentGenotypes[numParents] = latentParent->latentGenotype(locus)->value();
            allelesLostCounter_[numParents] += GeneticsImpl::falseNegativeCount(parentGenotypes[numParents], childGenotype);
            base = GeneticsImpl::shared(base, parentGenotypes[numParents]);// apply the latent parent genotype

            mutationFlag = mutationFlag.mutationMask(parentGenotypes[numParents]);
            if (mutationFlag.totalPositiveCount() > 0) {
//                fmt::print("Mutation flag is not zero: {}\n", core::io::serialize(mutationFlag));
                return -std::numeric_limits<Likelihood>::infinity();
            }

            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents + 1; ++parentIdx) {
                psIdxGen.reset(numParents + 1, parentIdx);

                // invert the latent parent immediately
                while (!psIdxGen.completed) {
                    base = childGenotype;

                    // Invert parent genotypes by the index of the selected subset
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // AND all (maybe inverted) parent genotypes plus the latent parent genotype
                    // with the inverted child genotype
                    for (size_t ii = 0; ii <= numParents; ++ii) {
                        base = GeneticsImpl::shared(base, parentGenotypes[ii]);
                    }

                    // Record the number of alleles that would be lost under that transmission
                    jointAllelesLostCounter_[jointIdx] += base.totalPositiveCount();
                    ++jointIdx;

                    // Invert again to get back to original state
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // Go to the next combination of parent indices
                    psIdxGen.next();
                }
            }
        }

        for (size_t ii = 0; ii < (unsigned long) std::pow(MAX_TRANSMISSIONS, numParents); ++ii) {
            auto kVec        = kVecs_[ii];// vector of numbers of transmission events for each parent
            kVec[numParents] = 1;         // the latent parent always has one transmission event

            double all_parents_lost_prob = 1.0;
            // Calculate the probability of losses in the transmission process
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = 1.0 - this->value()[kVec[parentIdx]];
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->value()(kVec[parentIdx]));
            }

            // losing the latentParent in a single transmission event
            all_parents_lost_prob *= 1.0 - this->value()[1];
            probs_[ii] += allelesLostCounter_[numParents] * std::log(1.0 - this->value()[1]);
            probs_[ii] += stp->value();

            // Calculate the probability of retentions in the transmission process where no alleles are lost
            probs_[ii] += jointAllelesLostCounter_[0] * std::log(1 - all_parents_lost_prob);
            size_t jointEventCounter = 1;

            for (size_t parentCount = 1; parentCount <= numParents; ++parentCount) {

                // Now we need to iterate over the subsets of events
                psIdxGen.reset(numParents + 1, parentCount);
                while (!psIdxGen.completed) {
                    // start with everything being multiplied out, i.e. all events happened
                    // for example (abc)
                    double tmpProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset
                        tmpProb /= (1.0 - this->value()[kVec[pIdx]]);
                    }
                    // what were left is the probability of the subset of parents retaining alleles
                    // i.e. (ab) if c was removed
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
//        if (std::isnan(llik)) {
//            llik = -std::numeric_limits<Likelihood>::infinity();
//        }


//        fmt::print("llik: {}\n", llik);
        if (std::isnan(llik) or llik <= -std::numeric_limits<Likelihood>::infinity()) {
//            fmt::print("NaN Encountered in SimpleLoss::calculateLogLikelihood with parent: {}\n", llik);
//            fmt::print("allelesLost: {}\n", core::io::serialize(allelesLostCounter_));
//            fmt::print("stp: {}\n", stp->value());
//            fmt::print("loss: {}\n", this->value()[1]);
//            fmt::print("retained: {}\n", 1.0 - this->value()[1]);
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

        return llik;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        auto loci                    = infection->loci();
        unsigned int allelesLost     = 0;
        unsigned int allelesRetained = 0;

        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype = infection->latentGenotype(locus)->value();
            allelesLost += GeneticsImpl::falseNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesRetained += GeneticsImpl::truePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);

            int anyMutation = GeneticsImpl::falsePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            if (anyMutation) {
                return -std::numeric_limits<Likelihood>::infinity();
            }
        }

        Likelihood allelesLostLlik = allelesLost * std::log(1.0 - this->value()[1]);
        Likelihood allelesRetainedLlik = allelesRetained * std::log(this->value()[1]);
        Likelihood stpLlik = stp->value();
        Likelihood llik = allelesLostLlik + allelesRetainedLlik + stpLlik;

//
//

//        fmt::print("llik: {}\n", llik);
        if (std::isnan(llik) or std::isinf(llik)) {
//            fmt::print("NaN Encountered in SimpleLoss::calculateLogLikelihood: {}\n", llik);
//            fmt::print("allelesLost: {}\n", allelesLost);
//            fmt::print("allelesRetained: {}\n", allelesRetained);
//            fmt::print("stp: {}\n", stp->value());
//            fmt::print("loss: {}\n", this->value()[1]);
//            fmt::print("retained: {}\n", 1.0 - this->value()[1]);
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

        return llik;
    }



    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE> parentGenotypes{};


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;

            auto mutationFlag = childGenotype;

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                mutationFlag = mutationFlag.mutationMask(parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            if (mutationFlag.totalPositiveCount() > 0) {
                return -std::numeric_limits<Likelihood>::infinity();
            }


            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents; ++parentIdx) {
                psIdxGen.reset(numParents, parentIdx);
                while (!psIdxGen.completed) {
                    base = childGenotype;

                    // Invert parent genotypes by the index of the selected subset
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // AND all (maybe inverted) parent genotypes with the inverted child genotype
                    for (size_t ii = 0; ii < numParents; ++ii) {
                        base = GeneticsImpl::shared(base, parentGenotypes[ii]);
                    }

                    // Record the number of alleles that would be lost under that transmission
                    jointAllelesLostCounter_[jointIdx] += base.totalPositiveCount();
                    ++jointIdx;

                    // Invert again to get back to original state
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // Go to the next combination of parent indices
                    psIdxGen.next();
                }
            }
        }

        for (size_t ii = 0; ii < (unsigned long) std::pow(MAX_TRANSMISSIONS, numParents); ++ii) {
            auto kVec = kVecs_[ii];

            double all_parents_lost_prob = 1.0;
            // Calculate the probability of losses in the transmission process
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = 1.0 - this->peek()[kVec[parentIdx]];
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->peek()(kVec[parentIdx]));
            }

            // Calculate the probability of retentions in the transmission process
            probs_[ii] += jointAllelesLostCounter_[0] * std::log(1 - all_parents_lost_prob);
            size_t jointEventCounter = 1;
            for (size_t parentCount = 1; parentCount < numParents; ++parentCount) {

                // Now we need to iterate over the subsets of events
                psIdxGen.reset(numParents, parentCount);
                while (!psIdxGen.completed) {
                    // start with everything being multiplied out, i.e. all events happened
                    // for example (abc)
                    double tmpProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset
                        tmpProb /= (1.0 - this->peek()[kVec[pIdx]]);
                    }
                    // what were left is the probability of the subset of parents retaining alleles
                    // i.e. (ab) if c was removed
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);

        if (std::isnan(llik)) {
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

        return llik;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE + 1> parentGenotypes{};



        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;

            auto mutationFlag = childGenotype;

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                mutationFlag = mutationFlag.mutationMask(parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            // calculate the latentParent contrib
            parentGenotypes[numParents] = latentParent->latentGenotype(locus)->value();
            allelesLostCounter_[numParents] += GeneticsImpl::falseNegativeCount(parentGenotypes[numParents], childGenotype);
            base = GeneticsImpl::shared(base, parentGenotypes[numParents]);// apply the latent parent genotype

            mutationFlag = mutationFlag.mutationMask(parentGenotypes[numParents]);
            if (mutationFlag.totalPositiveCount() > 0) {
                return -std::numeric_limits<Likelihood>::infinity();
            }

            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents + 1; ++parentIdx) {
                psIdxGen.reset(numParents + 1, parentIdx);

                // invert the latent parent immediately
                while (!psIdxGen.completed) {
                    base = childGenotype;

                    // Invert parent genotypes by the index of the selected subset
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // AND all (maybe inverted) parent genotypes plus the latent parent genotype
                    // with the inverted child genotype
                    for (size_t ii = 0; ii <= numParents; ++ii) {
                        base = GeneticsImpl::shared(base, parentGenotypes[ii]);
                    }

                    // Record the number of alleles that would be lost under that transmission
                    jointAllelesLostCounter_[jointIdx] += base.totalPositiveCount();
                    ++jointIdx;

                    // Invert again to get back to original state
                    for (const auto& kk : psIdxGen.curr) {
                        parentGenotypes[kk] = GeneticsImpl::invert(parentGenotypes[kk]);
                    }

                    // Go to the next combination of parent indices
                    psIdxGen.next();
                }
            }
        }

        for (size_t ii = 0; ii < (unsigned long) std::pow(MAX_TRANSMISSIONS, numParents); ++ii) {
            auto kVec        = kVecs_[ii];// vector of numbers of transmission events for each parent
            kVec[numParents] = 1;         // the latent parent always has one transmission event

            double all_parents_lost_prob = 1.0;
            // Calculate the probability of losses in the transmission process
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = 1.0 - this->peek()[kVec[parentIdx]];
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->peek()(kVec[parentIdx]));
            }

            // losing the latentParent in a single transmission event
            all_parents_lost_prob *= 1.0 - this->peek()[1];
            probs_[ii] += allelesLostCounter_[numParents] * std::log(1.0 - this->peek()[1]);
            probs_[ii] += stp->peek();

            // Calculate the probability of retentions in the transmission process where no alleles are lost
            probs_[ii] += jointAllelesLostCounter_[0] * std::log(1 - all_parents_lost_prob);
            size_t jointEventCounter = 1;

            for (size_t parentCount = 1; parentCount <= numParents; ++parentCount) {

                // Now we need to iterate over the subsets of events
                psIdxGen.reset(numParents + 1, parentCount);
                while (!psIdxGen.completed) {
                    // start with everything being multiplied out, i.e. all events happened
                    // for example (abc)
                    double tmpProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset
                        tmpProb /= (1.0 - this->peek()[kVec[pIdx]]);
                    }
                    // what were left is the probability of the subset of parents retaining alleles
                    // i.e. (ab) if c was removed
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
        if (std::isnan(llik)) {
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

        return llik;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        auto loci                    = infection->loci();
        unsigned int allelesLost     = 0;
        unsigned int allelesRetained = 0;

        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype = infection->latentGenotype(locus)->value();
            allelesLost += GeneticsImpl::falseNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesRetained += GeneticsImpl::truePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);

            int anyMutation = GeneticsImpl::falsePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            if (anyMutation) {
                return -std::numeric_limits<Likelihood>::infinity();
            }
        }

        Likelihood llik = allelesLost * std::log(1.0 - this->peek()[1]);
        llik += allelesRetained * std::log(this->peek()[1]);
        llik += stp->peek();

        if (std::isnan(llik)) {
            llik = -std::numeric_limits<Likelihood>::infinity();
        }

        return llik;
    }
;

}// namespace transmission_nets::model::transmission_process

#endif//TRANSMISSION_NETWORKS_APP_SIMPLELOSS_H
