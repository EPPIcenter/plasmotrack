//
// Created by Maxwell Murphy on 11/14/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_SEQUENTIALSCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SEQUENTIALSCHEDULER_H

#include "config.h"
namespace transmission_nets::impl::Model {

    template<typename T, typename Engine = boost::random::mt19937, typename Scheduler = core::samplers::Scheduler>
    struct SequentialSampleScheduler {
        SequentialSampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r);
        void step(const bool adapt = true) {
            scheduler_.step(adapt);
        }

        std::shared_ptr<State> state_{};
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> r_;
        Scheduler scheduler_;
    };

    template<typename T, typename Engine, typename Scheduler>
    SequentialSampleScheduler<T, Engine, Scheduler>::SequentialSampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r) : state_(std::move(state)), target_(std::move(target)), r_(r), scheduler_() {
        using namespace core::samplers;

        const double totalInfections = state_->infections.size();
        const double totalLoci = state_->loci.size();

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanCOI, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_coi",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 100});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanStrainsTransmitted, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_strains_tx",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 100});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->parentSetSizeProb, target_, 0, 1, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/parent_set_size_prob",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 100});

        int infection_idx_ = 0;
        for (const auto& infection : state_->infections) {

            const bool isSymptomatic = infection->isSymptomatic();
            const std::string infection_id = infection->id();
            const double upperBound = isSymptomatic ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size();

            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target_, 1.0, upperBound, r, 1, .1, 2),
                                        .id = fmt::format("ConstrainedContinuousRandomWalk/infection_duration/{}", infection_id),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = totalInfections * 10});
            if (!state_->null_model_) {
                scheduler_.registerSampler({//                    .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection_id], state_->latentParents[infection_idx_], target_, r),
                                            .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection_id], state_->latentParents[infection_idx_], infection->infectionDuration(), target_, r, 1.0, infection->isSymptomatic() ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size()),
                                            .id = fmt::format("JointGeneticsTimeSampler/infection/{}", infection_id),
                                            .weight = totalLoci});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler3<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection_id], state_->latentParents[infection_idx_], target_, r),
                                            .id = fmt::format("RandomAllelesBitSetSampler3/latent_genotypes/{}", infection_id),
                                            .weight = totalLoci});
                // Move the infection to another source in one step, and keep its latent parent able to be that source.
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::ParentRelinkSampler<T, Engine, InfectionEvent>>(infection, state_->latentParents[infection_idx_], state_->allowedRelationships->allowedParents(infection), target_, r, 1.0, upperBound, isSymptomatic ? state_->symptomaticInfectionDurationDist->value() : state_->asymptomaticInfectionDurationDist->value()),
                                            .id = fmt::format("ParentRelinkSampler/infection/{}", infection_id),
                                            .weight = totalLoci});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::LatentParentCoverSampler<T, Engine, InfectionEvent>>(infection, state_->latentParents[infection_idx_], target_, r, MAX_COI),
                                            .id = fmt::format("LatentParentCoverSampler/latent_parent_genotypes/{}", state_->latentParents[infection_idx_]->id()),
                                            .weight = totalLoci});
                // Move the infection to the other side of a clone in the ordering, and exchange their latent parents.
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>>(infection, state_->latentParents[infection_idx_], specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>::mutualPartners(state_->infections, state_->latentParents, *state_->allowedRelationships, infection_idx_), target_, r, 1.0, upperBound, isSymptomatic ? state_->symptomaticInfectionDurationDist->value() : state_->asymptomaticInfectionDurationDist->value()),
                                            .id = fmt::format("CloneOrderSwapSampler/ordering/{}", infection_id),
                                            .weight = totalLoci});
                for (const auto& [locus_label, locus] : state_->loci) {
                    if (infection->latentGenotype().contains(locus)) {
                        auto latentGenotype = infection->latentGenotype(locus);
                        // scheduler_.registerSampler({//                            .sampler = std::make_unique<genetics::RandomAllelesBitSetSampler2<T, Engine, GeneticsImpl, LocusImpl, ParentSetImpl>>(latentGenotype, locus, state_->parentSetList[infection_id], target_, r, MAX_COI),
                        // .sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r, MAX_COI),
                        // .id = fmt::format("Genotype {} {}", infection_id, locus->label),
                        // .weight = 5});
                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler4<T, Engine, InfectionEvent, GeneticsImpl, LocusImpl>>(infection, state_->latentParents[infection_idx_], locus, state_->allowedRelationships, target_, r, MAX_COI),
                                                    .id = fmt::format("RandomAllelesBitSetSampler4/latent_genotypes/{}/{}", infection_id, locus_label),
                                                    .weight = 5});
                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl, 1>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                        // scheduler_.registerSampler({.sampler = std::make_unique<genetics::SequentialAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                        //                             .weight  = 1});
                    }
                }
            }
            infection_idx_++;
        }


        if (!state_->null_model_) {
            for (const auto& infection : state_->latentParents) {
                for (const auto& [locus_label, locus] : state_->loci) {
                    if (infection->latentGenotype().contains(locus)) {
                        auto latentGenotype = infection->latentGenotype(locus);
                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r, MAX_COI),
                                                    .id = fmt::format("RandomAllelesBitSetSampler/latent_parent_genotypes/{}/{}", infection->id(), locus_label),
                                                    .weight = 5});

                        // scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl, 2>>(latentGenotype, target_, r),
                        //                             .weight = 1});
                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::SequentialAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                    }
                }
            }

            for (size_t i = 0; i < state_->expectedFalseNegatives.size(); ++i) {
                scheduler_.registerSampler({
                        .sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->expectedFalseNegatives[i], target_, 1e-6, .5, r, 1, .1, 2),
                        .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_negatives/{}", state_->infections[i]->id()),
                        .adaptationStart = 20,
                        .adaptationEnd = 200,
                        //                                        .update_start_     = 100,
                        .weight = 100
                        //                                        .weight          = 1
                });
            }

            for (size_t i = 0; i < state_->expectedFalsePositives.size(); ++i) {
                scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->expectedFalsePositives[i], target_, 1e-6, .5, r, 1, .1, 2),
                                            .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_positives/{}", state_->infections[i]->id()),
                                            .adaptationStart = 20,
                                            .adaptationEnd = 200,
                                            .weight = 100});
            }

            for (const auto& [locus_label, locus] : state_->loci) {
                auto alleleFrequencySampler = std::make_unique<SALTSampler<T, Engine>>(state_->alleleFrequencies->alleleFrequencies(locus), target_, r, 1, .1, 2);
                alleleFrequencySampler->setLowerLimit_(MIN_ALLELE_FREQUENCY);
                scheduler_.registerSampler({.sampler = std::move(alleleFrequencySampler),
                                            .id = fmt::format("SALTSampler/allele_frequencies/{}", locus_label),
                                            .adaptationStart = 20,
                                            .adaptationEnd = 200,
                                            .weight = 50});
            }
        }

    }
}// namespace transmission_nets::impl::ModelNine

#endif//TRANSMISSION_NETWORKS_APP_SEQUENTIALSCHEDULER_H
