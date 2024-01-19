//
// Created by Maxwell Murphy on 7/11/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLELOSSMUTATION_H
#define TRANSMISSION_NETWORKS_APP_SIMPLELOSSMUTATION_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Parameter.h"
#include "core/utils/numerics.h"
#include "model/transmission_process/NetworkBasedTransmissionProcess.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"

#include <memory>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    /*
     * Functional Node -- implements calculateLogLikelihood(child, parent_set)
     * Transmission process is a function of the loss and mutation rate per transmission event.
     * Internally calculates a 2x2 transition matrix M representing the probability of a gain or loss of a single allele, which
     * is then integrated over MAX_TRANSMISSIONS to get the probability of an allele being lost or gained
     */
    class SimpleLossMutation : public core::computation::Computation<std::array<double, (MAX_TRANSMISSIONS + 1) * 4>>,
                       public core::abstract::Observable<SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>,
                       public core::abstract::Cacheable<SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>,
                       public core::abstract::Checkpointable<SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>, std::array<double, (MAX_TRANSMISSIONS + 1) * 4>> {

        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;

    public:
        explicit SimpleLossMutation(p_Parameterdouble loss_prob, p_Parameterdouble mutation_rate, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb);

        std::array<double, (MAX_TRANSMISSIONS + 1) * 4> value() noexcept override;

        double getLossProb(int generation);
        double getMutationProb(int generation);
        double peekLossProb(int generation);
        double peekMutationProb(int generation);

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
        friend class core::abstract::Checkpointable<SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>, std::array<double, (MAX_TRANSMISSIONS + 1) * 4>>;
        friend class core::abstract::Cacheable<SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>>;


        std::array<std::array<unsigned int, MAX_PARENTSET_SIZE + 1>, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> kVecs_ = core::utils::initKvecs<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE>();
        p_Parameterdouble mutProb_;
        p_Parameterdouble lossProb_;
        std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb_;

        // Private store for probabilities when calculating the likelihood.
        std::array<double, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> probs_{};

        // tracks [x, y, z] where x, y, and z are the number of alleles that are lost in the transmission process.
        std::array<unsigned int, MAX_PARENTSET_SIZE + 1> allelesLostCounter_{0};

        // For example, parentset size = 3, tracks [a, b, c, ab, bc, ac, abc] where a, b, c, ab, bc, ac, and abc are the number of alleles that are lost in the transmission process.
        std::array<unsigned int, core::utils::const_pow(2, MAX_PARENTSET_SIZE + 1)> jointAllelesLostCounter_{0};
    };

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::SimpleLossMutation(p_Parameterdouble loss_prob, p_Parameterdouble mutation_rate, std::shared_ptr<InterTransmissionProbImpl> interTransmissionProb) : mutProb_(std::move(mutation_rate)), lossProb_(std::move(loss_prob)), interTransmissionProb_(std::move(interTransmissionProb)) {

        mutProb_->template registerCacheableCheckpointTarget(this);
        mutProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        lossProb_->template registerCacheableCheckpointTarget(this);
        lossProb_->add_post_change_listener([this]() {
            this->setDirty();
        });

        interTransmissionProb_->template registerCacheableCheckpointTarget(this);
        interTransmissionProb_->add_set_dirty_listener([this]() {
            this->setDirty();
        });

        this->value_.fill(0.0);

        this->setDirty();
        this->value();
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    std::array<double, (MAX_TRANSMISSIONS + 1) * 4> SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::value() noexcept {
        if (this->isDirty()) {
            this->value_.fill(0.0);

            /*
             * Sketch of transmission matrix and the indices
             *         0   1
             *       |-------
             *     0 | 0   2
             *     1 | 1   3
             *
             */

            // 0 -> 0
            this->value_[0] = 1.0 - this->mutProb_->value();
            // 1 -> 0
            this->value_[1] = this->lossProb_->value();
            // 0 -> 1
            this->value_[2] = this->mutProb_->value();
            // 1 -> 1
            this->value_[3] = 1.0 - this->lossProb_->value();

            for (unsigned int i = 1; i <= MAX_TRANSMISSIONS; i++) {
                unsigned int j = i - 1;
                // 0 -> 0
                this->value_[i * 4] = this->value_[j * 4] * this->value_[0] + this->value_[j * 4 + 2] * this->value_[1];
                // 1 -> 0
                this->value_[i * 4 + 1] = this->value_[j * 4 + 1] * this->value_[0] + this->value_[j * 4 + 3] * this->value_[1];
                // 0 -> 1
                this->value_[i * 4 + 2] = this->value_[j * 4] * this->value_[2] + this->value_[j * 4 + 2] * this->value_[3];
                // 1 -> 1
                this->value_[i * 4 + 3] = this->value_[j * 4 + 1] * this->value_[2] + this->value_[j * 4 + 3] * this->value_[3];
            }
            this->setClean();
        }
        return this->value_;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    double SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::getLossProb(int generation) {
        return this->value()[(generation - 1) * 4 + 1];
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    double SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::getMutationProb(int generation) {
        return this->value()[(generation - 1) * 4 + 2];
    }


    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    double SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLossProb(int generation) {
        return this->peek()[(generation - 1) * 4 + 1];
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    double SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekMutationProb(int generation) {
        return this->peek()[(generation - 1) * 4 + 2];
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE> parentGenotypes{};
        int totalNegativeAlleles = 0;


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;
            totalNegativeAlleles += childGenotype.totalNegativeCount();

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                ++parentIdx;
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
            // Calculate the probability of losses in the transmission process -- i.e. 1->0 or 0->0
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = this->getLossProb(kVec[parentIdx]);
                const double not_mutated_prob = 1.0 - this->getMutationProb(kVec[parentIdx]);
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob) +
                              (totalNegativeAlleles - allelesLostCounter_[parentIdx]) * std::log(not_mutated_prob);

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
                    double eventProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset and replace with a 0->0 event
                        eventProb /= this->getLossProb(kVec[pIdx]);
                        eventProb *= (1.0 - this->getMutationProb(kVec[pIdx]));
                    }
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - eventProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        return core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE + 1> parentGenotypes{};
        int totalNegativeAlleles = 0;


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;
            totalNegativeAlleles += childGenotype.totalNegativeCount();

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            // calculate the latentParent contrib
            parentGenotypes[numParents] = latentParent->latentGenotype(locus)->value();
            allelesLostCounter_[numParents] += GeneticsImpl::falseNegativeCount(parentGenotypes[numParents], childGenotype);
            base = GeneticsImpl::shared(base, parentGenotypes[numParents]); // apply the latent parent genotype

//            // only allow loss from the source parent, must be a complete subset of the child genotype
//            auto mutationFlag = childGenotype.mutationMask(parentGenotypes[numParents]);
//            if (mutationFlag.totalPositiveCount() > 0) {
//                return -std::numeric_limits<Likelihood>::infinity();
//            }

            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents + 1; ++parentIdx) {
                psIdxGen.reset(numParents + 1, parentIdx);

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
                const double loss_prob = this->getLossProb(kVec[parentIdx]);
                const double not_mutated_prob = 1.0 - this->getMutationProb(kVec[parentIdx]);
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob) +
                              (totalNegativeAlleles - allelesLostCounter_[parentIdx]) * std::log(not_mutated_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->value()(kVec[parentIdx]));
            }

            // losing the latentParent in a single transmission event
            const double lost_in_one_gen = this->getLossProb(1);
            const double not_mutated_in_one_gen = 1.0 - this->getMutationProb(1);
            all_parents_lost_prob *= lost_in_one_gen;
            probs_[ii] += allelesLostCounter_[numParents] * std::log(lost_in_one_gen) +
                          (totalNegativeAlleles - allelesLostCounter_[numParents]) * std::log(not_mutated_in_one_gen);
            // tack on probability of latentParent being drawn from the background population
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
                        tmpProb /= this->getLossProb(kVec[pIdx]);
                        tmpProb *= (1.0 - this->getMutationProb(kVec[pIdx]));
                    }
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
        return llik;
    }


    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::calculateLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        auto loci                    = infection->loci();
        unsigned int allelesLost     = 0;
        unsigned int allelesRetained = 0;
        unsigned int allelesMutated  = 0;
        unsigned int allelesNotMutated = 0;

        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype = infection->latentGenotype(locus)->value();
            allelesLost += GeneticsImpl::falseNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesRetained += GeneticsImpl::truePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesMutated += GeneticsImpl::falsePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesNotMutated += GeneticsImpl::trueNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
//            if (allelesMutated > 0) {
//                return -std::numeric_limits<Likelihood>::infinity();
//            }
        }

        // 1 -> 0
        Likelihood llik = allelesLost * std::log(this->getLossProb(1));
        // 1 -> 1
        llik += allelesRetained * std::log(1 - this->getLossProb(1));
        // 0 -> 1
        llik += allelesMutated * std::log(this->getMutationProb(1));
        // 0 -> 0
        llik += allelesNotMutated * std::log(1 - this->getMutationProb(1));
        // Add in the probability of the transmission process
        llik += stp->value();

        return llik;
    }



    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE> parentGenotypes{};
        int totalNegativeAlleles = 0;


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;
            totalNegativeAlleles += childGenotype.totalNegativeCount();

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                ++parentIdx;
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
            // Calculate the probability of losses in the transmission process -- i.e. 1->0 or 0->0
            for (size_t parentIdx = 0; parentIdx < numParents; ++parentIdx) {
                // probability that an allele is lost after k generations
                const double loss_prob = this->peekLossProb(kVec[parentIdx]);
                const double not_mutated_prob = 1.0 - this->peekMutationProb(kVec[parentIdx]);
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob) +
                              (totalNegativeAlleles - allelesLostCounter_[parentIdx]) * std::log(not_mutated_prob);

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
                    double eventProb = all_parents_lost_prob;
                    for (const auto& pIdx : psIdxGen.curr) {
                        // remove each parent inside the indexed subset and replace with a 0->0 event
                        eventProb /= this->peekLossProb(kVec[pIdx]);
                        eventProb *= (1.0 - this->peekMutationProb(kVec[pIdx]));
                    }
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - eventProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        return core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, const core::containers::ParentSet<core::containers::Infection<GeneticsImpl>>& parentSet, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        size_t numParents = parentSet.size();
        auto loci         = infection->loci();

        core::utils::generators::CombinationIndicesGenerator psIdxGen;
        // re-initialize the arrays
        probs_.fill(0);
        allelesLostCounter_.fill(0);
        jointAllelesLostCounter_.fill(0);
        std::array<GeneticsImpl, MAX_PARENTSET_SIZE + 1> parentGenotypes{};
        int totalNegativeAlleles = 0;


        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype         = infection->latentGenotype(locus)->value();
            auto invertedChildGenotype = GeneticsImpl::invert(childGenotype);
            auto base                  = invertedChildGenotype;
            size_t parentIdx           = 0;
            totalNegativeAlleles += childGenotype.totalNegativeCount();

            // iterate over all parents to get their genotypes and calculate the number of
            // alleles that are lost in the transmission process.
            for (const auto& parent : parentSet) {
                parentGenotypes[parentIdx] = parent->latentGenotype(locus)->value();
                allelesLostCounter_[parentIdx] += GeneticsImpl::falseNegativeCount(parentGenotypes[parentIdx], childGenotype);
                base = GeneticsImpl::shared(base, parentGenotypes[parentIdx]);
                ++parentIdx;
            }

            // calculate the latentParent contrib
            parentGenotypes[numParents] = latentParent->latentGenotype(locus)->value();
            allelesLostCounter_[numParents] += GeneticsImpl::falseNegativeCount(parentGenotypes[numParents], childGenotype);
            base = GeneticsImpl::shared(base, parentGenotypes[numParents]); // apply the latent parent genotype

//            // only allow loss from the source parent, must be a complete subset of the child genotype
//            auto mutationFlag = childGenotype.mutationMask(parentGenotypes[numParents]);
//            if (mutationFlag.totalPositiveCount() > 0) {
//                return -std::numeric_limits<Likelihood>::infinity();
//            }

            // Calculating the alternative situation, where alleles are not lost
            jointAllelesLostCounter_[0] += base.totalPositiveCount();
            size_t jointIdx = 1;
            for (parentIdx = 1; parentIdx <= numParents + 1; ++parentIdx) {
                psIdxGen.reset(numParents + 1, parentIdx);

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
                const double loss_prob = this->peekLossProb(kVec[parentIdx]);
                const double not_mutated_prob = 1.0 - this->peekMutationProb(kVec[parentIdx]);
                all_parents_lost_prob *= loss_prob;
                probs_[ii] += allelesLostCounter_[parentIdx] * std::log(loss_prob) +
                              (totalNegativeAlleles - allelesLostCounter_[parentIdx]) * std::log(not_mutated_prob);

                // tack on probability of k generations for parent at parentIdx
                probs_[ii] += std::log(interTransmissionProb_->peek()(kVec[parentIdx]));
            }

            // losing the latentParent in a single transmission event
            const double lost_in_one_gen = this->peekLossProb(1);
            const double not_mutated_in_one_gen = 1.0 - this->peekMutationProb(1);
            all_parents_lost_prob *= lost_in_one_gen;
            probs_[ii] += allelesLostCounter_[numParents] * std::log(lost_in_one_gen) +
                          (totalNegativeAlleles - allelesLostCounter_[numParents]) * std::log(not_mutated_in_one_gen);
            // tack on probability of latentParent being drawn from the background population
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
                        tmpProb /= this->peekLossProb(kVec[pIdx]);
                        tmpProb *= (1.0 - this->peekMutationProb(kVec[pIdx]));
                    }
                    probs_[ii] += jointAllelesLostCounter_[jointEventCounter] * std::log(1 - tmpProb);
                    ++jointEventCounter;
                    psIdxGen.next();
                }
            }
        }
        Likelihood llik = core::utils::logSumExp(probs_.begin(), probs_.begin() + (unsigned int) std::pow(MAX_TRANSMISSIONS, numParents) - 1);
        return llik;
    }


    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE, typename InterTransmissionProbImpl, typename SourceTransmissionProcessImpl>
    template<typename GeneticsImpl>
    Likelihood SimpleLossMutation<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE, InterTransmissionProbImpl, SourceTransmissionProcessImpl>::peekLogLikelihood(std::shared_ptr<core::containers::Infection<GeneticsImpl>> infection, std::shared_ptr<core::containers::Infection<GeneticsImpl>> latentParent, std::shared_ptr<SourceTransmissionProcessImpl> stp) {
        auto loci                    = infection->loci();
        unsigned int allelesLost     = 0;
        unsigned int allelesRetained = 0;
        unsigned int allelesMutated  = 0;
        unsigned int allelesNotMutated = 0;

        // Calculate the number of events that occur in the transmission process, where
        // events are either the loss or retention of alleles.
        for (const auto& locus : loci) {
            auto childGenotype = infection->latentGenotype(locus)->value();
            allelesLost += GeneticsImpl::falseNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesRetained += GeneticsImpl::truePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesMutated += GeneticsImpl::falsePositiveCount(latentParent->latentGenotype(locus)->value(), childGenotype);
            allelesNotMutated += GeneticsImpl::trueNegativeCount(latentParent->latentGenotype(locus)->value(), childGenotype);
//            if (allelesMutated > 0) {
//                return -std::numeric_limits<Likelihood>::infinity();
//            }
        }

        // 1 -> 0
        Likelihood llik = allelesLost * std::log(this->peekLossProb(1));
        // 1 -> 1
        llik += allelesRetained * std::log(1 - this->peekLossProb(1));
        // 0 -> 1
        llik += allelesMutated * std::log(this->peekMutationProb(1));
        // 0 -> 0
        llik += allelesNotMutated * std::log(1 - this->peekMutationProb(1));

        // Add in the probability of the transmission process
        llik += stp->peek();

        return llik;
    }



}



#endif//TRANSMISSION_NETWORKS_APP_SIMPLELOSSMUTATION_H
