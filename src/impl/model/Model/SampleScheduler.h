//
// Created by mmurphy on 10/29/21.
//
#ifndef TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H

#include "config.h"

#include "core/samplers/scheduler/RandomizedScheduler.h"
#include "core/samplers/genetics/AlleleSwapSampler.h"
#include "core/samplers/genetics/StrainCountSampler.h"
#include "core/samplers/specialized/CloneOrderSwapSampler.h"
#include "core/samplers/specialized/InformedOrderJumpSampler.h"
#include "core/samplers/specialized/JointGeneticsTimeSampler.h"
#include "core/samplers/specialized/ParentRelinkSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler4.h"
#include "core/samplers/general/SALTSampler.h"


#include <core/samplers/general/ConstrainedContinuousRandomWalk.h>
namespace transmission_nets::impl::Model {

    template<typename T, typename Engine = boost::random::mt19937, typename Scheduler = core::samplers::RandomizedScheduler<Engine>>
    struct SampleScheduler {
        // Without `updatesPerIteration`, each Sampling iteration makes one update per registered sampler; that is also
        // how replica exchange builds it.
        SampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r, std::optional<int> updatesPerIteration = std::nullopt);
        void step(const bool adapt = true) {
            scheduler_.step(adapt);
        }

        std::shared_ptr<State> state_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> r_;
        Scheduler scheduler_;
    };

    template<typename T, typename Engine, typename Scheduler>
    SampleScheduler<T, Engine, Scheduler>::SampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r, std::optional<int> updatesPerIteration) : state_(std::move(state)), target_(std::move(target)), r_(r), scheduler_(r_, updatesPerIteration) {
        using namespace core::samplers;

        [[maybe_unused]] double totalInfections = state_->infections.size();
        [[maybe_unused]] double totalLoci = state_->loci.size();

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanCOI, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_coi",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanStrainsTransmitted, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_strains_tx",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->parentSetSizeProb, target_, 0, 1, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/parent_set_size_prob",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5});

        int infection_idx_ = 0;
        for (auto& infection : state_->infections) {
            const bool isSymptomatic = infection->isSymptomatic();
            double upperBound = isSymptomatic ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size();

            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target_, 1.0, upperBound, r, 1, .1, 2),
                                        .id = fmt::format("ConstrainedContinuousRandomWalk/infection_duration/{}", infection->id()),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 5});
            if (!state_->null_model_) {
                scheduler_.registerSampler({//                    .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], target_, r),
                                            .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], infection->infectionDuration(), target_, r, 1.0, infection->isSymptomatic() ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size()),
                                            .id = fmt::format("JointGeneticsTimeSampler/infection/{}", infection->id()),
                                            .weight = 5});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler3<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], target_, r),
                                            .id = fmt::format("RandomAllelesBitSetSampler3/latent_genotypes/{}", infection->id()),
                                            .weight = 5});
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::ParentRelinkSampler<T, Engine, InfectionEvent>>(infection, state_->latentParents[infection_idx_], state_->allowedRelationships->allowedParents(infection), target_, r, 1.0, upperBound, isSymptomatic ? state_->symptomaticInfectionDurationDist->value() : state_->asymptomaticInfectionDurationDist->value(), 0.01, 0.01, RELINK_TRIES),
                                            .id = fmt::format("ParentRelinkSampler/infection/{}", infection->id()),
                                            .weight = 5});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::StrainCountSampler<T, Engine, InfectionEvent>>(infection, target_, r, MAX_COI),
                                            .id = fmt::format("StrainCountSampler/latent_genotypes/{}", infection->id()),
                                            .weight = 5});
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>>(infection, specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>::mutualPartners(state_->infections, *state_->allowedRelationships, infection_idx_), target_, r, 1.0, upperBound, isSymptomatic ? state_->symptomaticInfectionDurationDist->value() : state_->asymptomaticInfectionDurationDist->value()),
                                            .id = fmt::format("CloneOrderSwapSampler/ordering/{}", infection->id()),
                                            .weight = 5});
                {
                    std::vector<std::shared_ptr<InfectionEvent>> partners;
                    for (const auto& partner : specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>::mutualPartners(state_->infections, *state_->allowedRelationships, infection_idx_)) {
                        partners.push_back(partner.infection);
                    }
                    scheduler_.registerSampler({.sampler = std::make_unique<specialized::InformedOrderJumpSampler<T, Engine, InfectionEvent>>(infection, std::move(partners), target_, r, 1.0, upperBound, isSymptomatic ? state_->symptomaticInfectionDurationDist->value() : state_->asymptomaticInfectionDurationDist->value(), INFORMED_JUMP_MAX_TRIES),
                                                .id = fmt::format("InformedOrderJumpSampler/ordering/{}", infection->id()),
                                                .weight = 5});
                }
                for (const auto& [locus_label, locus] : state_->loci) {
                    if (infection->latentGenotype().contains(locus)) {
                        auto latentGenotype = infection->latentGenotype(locus);
                        //                     scheduler_.registerSampler({
                        // //                            .sampler = std::make_unique<genetics::RandomAllelesBitSetSampler2<T, Engine, GeneticsImpl, LocusImpl, ParentSetImpl>>(latentGenotype, locus, state_->parentSetList[infection->id()], target_, r, MAX_COI),
                        //                             .sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r, MAX_COI),
                        //                             .id      = fmt::format("Genotype {} {}", infection->id(), locus->label),
                        //                             .weight = 5
                        //                     });
                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler4<T, Engine, InfectionEvent, GeneticsImpl, LocusImpl>>(infection, state_->latentParents[infection_idx_], locus, state_->allowedRelationships, target_, r, MAX_COI),
                                                    .id = fmt::format("RandomAllelesBitSetSampler4/latent_genotypes/{}/{}", infection->id(), locus->label),
                                                    .weight = 5});
                        // Swaps pay where the observation says little about which allele is present: chosen as often as other
                        // per-locus samplers at loci with no allele observed, a tenth as often at loci with several, and never at loci with one.
                        const unsigned int observedAlleles = infection->observedGenotype().contains(locus) ? infection->observedGenotype(locus)->value().totalPositiveCount() : 0;
                        if (observedAlleles != 1) {
                            scheduler_.registerSampler({.sampler = std::make_unique<genetics::AlleleSwapSampler<T, Engine, GeneticsImpl>>(infection->latentGenotype(locus), target_, r),
                                                    .id = fmt::format("AlleleSwapSampler/latent_genotypes/{}/{}", infection->id(), locus->label),
                                                    .weight = observedAlleles == 0 ? 5.0 : 0.5});
                        }

                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl, 1>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::SequentialAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                    }
                }
                infection_idx_++;
            }
        }





        for (size_t i = 0; i < state_->expectedFalseNegatives.size(); ++i) {
            scheduler_.registerSampler({
                    .sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->expectedFalseNegatives[i], target_, 1e-6, .5, r, 1, .1, 2),
                    .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_negatives/{}", state_->infections[i]->id()),
                    .adaptationStart = 20,
                    .adaptationEnd = 200,
                    //                                        .update_start_     = 100,
                    .weight = 5
                    //                                        .weight          = 1
            });
        }

        for (size_t i = 0; i < state_->expectedFalsePositives.size(); ++i) {
            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->expectedFalsePositives[i], target_, 1e-6, .5, r, 1, .1, 2),
                                        .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_positives/{}", state_->infections[i]->id()),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 5});
        }

        for ([[maybe_unused]] const auto& [locus_label, locus] : state_->loci) {
            auto alleleFrequencySampler = std::make_unique<SALTSampler<T, Engine>>(state_->alleleFrequencies->alleleFrequencies(locus), target_, r, 1, .1, 2);
            alleleFrequencySampler->setLowerLimit_(MIN_ALLELE_FREQUENCY);
            scheduler_.registerSampler({.sampler = std::move(alleleFrequencySampler),
                                        .id = fmt::format("SALTSampler/allele_frequencies/{}", locus->label),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 1});
        }
    }
}// namespace transmission_nets::impl::ModelNine

#endif//TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
