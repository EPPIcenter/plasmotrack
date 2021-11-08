//
// Created by mmurphy on 10/29/21.
//

#include "Model.h"

#include <utility>

namespace transmission_nets::impl::ModelFive {

    Model::Model(std::shared_ptr<State> state) : state_(std::move(state)) {
        likelihood.add_set_dirty_listener([=, this]() {
          this->setDirty();
        });
        likelihood.registerCacheableCheckpointTarget(this);

        intp = std::make_shared<InterTransmissionProbImpl>(state_->geometricGenerationProb);
        nodeTransmissionProcess = std::make_shared<NodeTransmissionImpl>(state_->mutationProb, state_->lossProb, intp);
        coiProb = std::make_shared<COIProbabilityImpl>(state_->meanCOI);

        // Register Priors
        likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->mutationProb, state_->mutationProbPriorAlpha, state_->mutationProbPriorBeta));
        likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->lossProb, state_->lossProbPriorAlpha, state_->lossProbPriorBeta));
        likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->mutationProb, state_->mutationProbPriorAlpha, state_->mutationProbPriorBeta));
        likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanCOI, state_->meanCOIPriorShape, state_->meanCOIPriorScale));
        likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->geometricGenerationProb, state_->geometricGenerationProbPriorAlpha, state_->geometricGenerationProbPriorBeta));
        //        likelihood.addTarget(new core::distributions::GammaLogPDF(state_->infectionDurationShape, state_->infectionDurationShapePriorShape, state_->infectionDurationShapePriorScale));
        //        likelihood.addTarget(new core::distributions::GammaLogPDF(state_->infectionDurationScale, state_->infectionDurationScalePriorShape, state_->infectionDurationScalePriorScale));
        for (auto& obs : state_->observationFalsePositiveRates) {
            likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(obs, state_->obsFPRPriorAlpha, state_->obsFPRPriorBeta));
        }
        for (auto& obs : state_->observationFalseNegativeRates) {
            likelihood.addTarget(std::make_shared<core::distributions::BetaLogPDF>(obs, state_->obsFNRPriorAlpha, state_->obsFNRPriorBeta));
        }

        int i = 0;
        for (auto &infection : state_->infections) {
            likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(infection->infectionDuration(), state_->infectionDurationShape, state_->infectionDurationScale));
            alleleCountAccumulators.push_back(std::make_shared<AlleleCounterAccumulator>());

            for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
                alleleCounters.push_back(std::make_shared<AlleleCounterImpl>(infection->latentGenotype(locus), obsGenotype));
                alleleCountAccumulators.back()->addTarget(alleleCounters.back());
            }

            observationProcessLikelihood = std::make_shared<model::observation_process::ObservationProcessLikelihood<AlleleCounterAccumulator>>(
                    alleleCountAccumulators.back(),
                    state_->observationFalseNegativeRates[i],
                    state_->observationFalsePositiveRates[i]);
            likelihood.addTarget(observationProcessLikelihood);
            parentSetList.push_back(std::make_shared<ParentSetImpl>(state_->infectionEventOrdering, infection, state_->allowedParents[infection]));
            i++;
        }

        for (unsigned int j = 0; j < state_->infections.size(); ++j) {
            auto infection = state_->infections[j];
            auto parentSet = parentSetList[j];

            sourceTransmissionProcessList.push_back(std::make_shared<SourceTransmissionImpl>(
                    coiProb,
                    state_->alleleFrequencies,
                    infection));

            transmissionProcessList.push_back(std::make_shared<TransmissionProcess>(
                    nodeTransmissionProcess,
                    sourceTransmissionProcessList.back(),
                    infection,
                    parentSet));

            likelihood.addTarget(transmissionProcessList.back());
        }
        this->setDirty();
    }

    Model::Model(State &state) : Model(std::make_shared<State>(state)){}

    std::string Model::identifier() {
        return "Model5";
    }

    Likelihood Model::value() {
        if (isDirty()) {
            value_ = likelihood.value();
            this->setClean();
        }
        return value_;
    }
};
