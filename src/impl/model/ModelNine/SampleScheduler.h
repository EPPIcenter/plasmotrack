//
// Created by mmurphy on 10/29/21.
//
#ifndef TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H

#include "config.h"

#include "core/samplers/scheduler/RandomizedScheduler.h"
#include "core/samplers/specialized/JointGeneticsTimeSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler4.h"
#include "core/samplers/general/SALTSampler.h"


#include <core/samplers/general/ConstrainedContinuousRandomWalk.h>
namespace transmission_nets::impl::ModelNine {

    template<typename T, typename Engine = boost::random::mt19937, typename Scheduler = core::samplers::RandomizedScheduler<Engine>>
    struct SampleScheduler {
        SampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r, int samplesPerStep);
        void step() {
            scheduler_.step();
        }

        std::shared_ptr<State> state_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> r_;
        Scheduler scheduler_;
    };

    template<typename T, typename Engine, typename Scheduler>
    SampleScheduler<T, Engine, Scheduler>::SampleScheduler(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r, int samplesPerStep) : state_(std::move(state)), target_(std::move(target)), r_(r), scheduler_(r_, samplesPerStep) {
        using namespace core::samplers;

        [[maybe_unused]] double totalInfections = state_->infections.size();
        [[maybe_unused]] double totalLoci = state_->loci.size();

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanCOI, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "Mean COI",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5,
                                    .debug = false});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanStrainsTransmitted, target_, 1.0, 20, r, 1, .1, 1),
                                    .id = "Mean Strains Tx",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5,
                                    .debug = false});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->parentSetSizeProb, target_, 0, 1, r, 1, .1, 1),
                                    .id = "Parent Set Size",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200,
                                    .weight = 5,
                                    .debug = false});

        int infection_idx_ = 0;
        for (auto& infection : state_->infections) {
            bool isSymptomatic = infection->isSymptomatic();
            double upperBound = isSymptomatic ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size();

            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target_, 1.0, upperBound, r, 1, .1, 2),
                                        .id = fmt::format("Infection Duration {}", infection->id()),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 5,
                                        .debug = false});
            if (!state_->null_model_) {
                scheduler_.registerSampler({//                    .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], target_, r),
                                            .sampler = std::make_unique<specialized::JointGeneticsTimeSampler<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], infection->infectionDuration(), target_, r, 1.0, infection->isSymptomatic() ? state_->symptomaticInfectionDurationDist->value().size() : state_->asymptomaticInfectionDurationDist->value().size()),
                                            .id = fmt::format("Infection Alleles/Infection Duration {}", infection->id()),
                                            .weight = 5});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler3<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, state_->parentSetList[infection->id()], state_->latentParents[infection_idx_], target_, r),
                                            .id = fmt::format("Infection Alleles {}", infection->id()),
                                            .weight = 5,
                                            .debug = false});
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
                                                    .id = fmt::format("Genotype4 {} {}", infection->id(), locus->label),
                                                    .weight = 5,
                                                    .debug = true});

                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl, 1>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                        //                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::SequentialAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                        //                                                                    .weight  = 1});
                    }
                }
                infection_idx_++;
            }
        }


        for (auto& infection : state_->latentParents) {
            for (const auto& [locus_label, locus] : state_->loci) {
                if (infection->latentGenotype().contains(locus)) {
                    auto latentGenotype = infection->latentGenotype(locus);
                    scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r, MAX_COI),
                                                .id = fmt::format("Latent Genotype {} {}", infection->id(), locus_label),
                                                .weight = 5,
                                                .debug = false});

                    // scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl, 2>>(latentGenotype, target_, r),
                    //                             .weight  = 1});
                    // scheduler_.registerSampler({.sampler = std::make_unique<genetics::SequentialAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                    //                             .weight  = 1, .debug = false});
                }
            }
        }



        for ([[maybe_unused]] auto& infFNR : state_->expectedFalseNegatives) {
            scheduler_.registerSampler({
                    .sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infFNR, target_, 1e-6, .5, r, 1, .1, 2),
                    .id = fmt::format("False Negative Rate"),
                    .adaptationStart = 20,
                    .adaptationEnd = 200,
                    //                                        .update_start_     = 100,
                    .weight = 5
                    //                                        .weight          = 1
            });
        }

        for ([[maybe_unused]] auto& infFPR : state_->expectedFalsePositives) {
            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infFPR, target_, 1e-6, .5, r, 1, .1, 2),
                                        .id = fmt::format("False Positive Rate"),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 5});
        }

        for ([[maybe_unused]] const auto& [locus_label, locus] : state_->loci) {
            scheduler_.registerSampler({.sampler = std::make_unique<SALTSampler<T, Engine>>(state_->alleleFrequencies->alleleFrequencies(locus), target_, r, 1, .1, 2),
                                        .id = fmt::format("Allele Freq {}", locus->label),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200,
                                        .weight = 1});
        }
    }
}// namespace transmission_nets::impl::ModelNine

#endif//TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
