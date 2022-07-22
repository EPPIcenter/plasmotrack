//
// Created by mmurphy on 10/29/21.
//

#include "config.h"

#ifndef TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H

namespace transmission_nets::impl::ModelSeven {

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

        long double totalInfections = state_->infections.size();
        long double totalLoci       = state_->loci.size();

        scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->geometricGenerationProb, target_, 0.0, 1.0, r, .01, .01, 2),
                                    .id              = "Geometric Generation Prob",
                                    .adaptationStart = 0,
                                    .adaptationEnd   = 2000,
                                    .weight          = totalInfections * 10,
                                    .debug           = true});

        scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->lossProb, target_, 0.0, 1.0, r, .01, .01, 2),
                                    .id              = "Loss Probability",
                                    .adaptationStart = 0,
                                    .adaptationEnd   = 2000,
                                    .weight          = totalInfections * 10,
                                    .debug           = true});

        scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state_->meanCOI, target_, 0.0, 100, r, .01, .01, 1),
                                    .id              = "Mean COI",
                                    .adaptationStart = 0,
                                    .adaptationEnd   = 2000,
                                    .weight          = totalInfections * 10,
                                    .debug           = true});

        for (auto& infection : state_->infections) {
            scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target_, 1.0, 1000.0, r, 1, .1, 100),
                                        .id              = fmt::format("Infection Duration {}", infection->id()),
                                        .adaptationStart = 20,
                                        .adaptationEnd   = 2000,
                                        .weight          = totalLoci * 100});
            for (const auto& [locus_label, locus] : state_->loci) {
                if (infection->latentGenotype().contains(locus)) {
                    auto latentGenotype = infection->latentGenotype(locus);
                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                                                                    .id      = fmt::format("Genotype {} {}", infection->id(), locus->label),
                                                                    //                                                .updateStart = 1000,
                                                                    .weight = 5});
//                    scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
//                                                .weight  = 1});
                }
            }
        }

        for (auto& infection : state_->latentParents) {
            for (const auto& [locus_label, locus] : state_->loci) {
                if (infection->latentGenotype().contains(locus)) {
                    auto latentGenotype = infection->latentGenotype(locus);
                                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
                                                                    .id      = fmt::format("Latent Genotype {} {}", infection->id(), locus->label),
                                                                    //                                                .updateStart = 1000,
                                                                    .weight = 5});
//                    scheduler_.registerSampler({.sampler = std::make_unique<genetics::ZanellaAllelesBitSetSampler<T, Engine, GeneticsImpl>>(latentGenotype, target_, r),
//                                                .weight  = 1});
                }
            }
        }

        for (auto& infFNR : state_->expectedFalseNegatives) {
            scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infFNR, target_, 0.0, .5, r, .1),
                                        .id              = fmt::format("False Negative Rate"),
                                        .adaptationStart = 20,
                                        .adaptationEnd   = 2000,
                                        .weight          = totalLoci * 10});
        }

        for (auto& infFPR : state_->expectedFalsePositives) {
            scheduler_.registerSampler({.sampler         = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infFPR, target_, 0.0, .5, r, .1),
                                        .id              = fmt::format("False Positive Rate"),
                                        .adaptationStart = 20,
                                        .adaptationEnd   = 2000,
                                        .weight          = totalLoci * 10});
        }

        for (const auto& [locus_label, locus] : state_->loci) {
            scheduler_.registerSampler({.sampler         = std::make_unique<SALTSampler<T, Engine>>(state_->alleleFrequencies->alleleFrequencies(locus), target_, r, 1, .01, 10),
                                        .id              = fmt::format("Allele Freq {}", locus->label),
                                        .adaptationStart = 20,
                                        .adaptationEnd   = 2000,
                                        .weight          = totalInfections});
        }
    }
}// namespace transmission_nets::impl::ModelSeven

#endif//TRANSMISSION_NETWORKS_APP_SAMPLESCHEDULER_H
