//
// Created by mmurphy on 10/29/21.
//

#include "Model.h"

#include <utility>

namespace transmission_nets::impl::ModelEight {

    Model::Model(std::shared_ptr<State> state, double temperature) : state_(std::move(state)), temperature(temperature) {
        likelihood.add_set_dirty_listener([=, this]() {
            this->setDirty();
        });
        likelihood.registerCacheableCheckpointTarget(this);

        intp                    = std::make_shared<InterTransmissionProbImpl>(state_->geometricGenerationProb);
        nodeTransmissionProcess = std::make_shared<NodeTransmissionImpl>(state_->lossProb, intp);
        coiProb                 = std::make_shared<COIProbabilityImpl>(state_->meanCOI);

        // Register Priors
        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->lossProb, state_->lossProbPriorAlpha, state_->lossProbPriorBeta));
        prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanCOI, state_->meanCOIPriorShape, state_->meanCOIPriorScale));
        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->geometricGenerationProb, state_->geometricGenerationProbPriorAlpha, state_->geometricGenerationProbPriorBeta));
        //        likelihood.addTarget(new core::distributions::GammaLogPDF(state_->infectionDurationShape, state_->infectionDurationShapePriorShape, state_->infectionDurationShapePriorScale));
        //        likelihood.addTarget(new core::distributions::GammaLogPDF(state_->infectionDurationScale, state_->infectionDurationScalePriorShape, state_->infectionDurationScalePriorScale));
        for (auto& obs : state_->expectedFalsePositives) {
            likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFPRPriorShape, state_->obsFPRPriorScale));
        }
        for (auto& obs : state_->expectedFalseNegatives) {
            likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFNRPriorShape, state_->obsFNRPriorScale));
        }

        int i = 0;
        for (auto& infection : state_->infections) {
            // infection duration likelihood

            if (infection->isSymptomatic()) {
                likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(infection->infectionDuration(), state_->symptomaticInfectionDurationShape, state_->symptomaticInfectionDurationScale));
            } else {
                likelihood.addTarget(std::make_shared<core::distributions::GammaLogPDF>(infection->infectionDuration(), state_->asymptomaticInfectionDurationShape, state_->asymptomaticInfectionDurationScale));
            }


            for (auto& [locus, obsGenotype] : infection->observedGenotype()) {
                observationProcessLikelihoodList.push_back(std::make_shared<ObservationProcessImpl>(
                        obsGenotype,
                        infection->latentGenotype(locus),
                        state_->expectedFalsePositives[i],
                        state_->expectedFalseNegatives[i]));
                likelihood.addTarget(observationProcessLikelihoodList.back());
            }

            state_->parentSetList[infection->id()] = std::make_shared<ParentSetImpl>(state_->infectionEventOrdering, infection, state_->allowedParents[infection]);
            i++;
        }

        for (unsigned int j = 0; j < state_->infections.size(); ++j) {
            auto infection = state_->infections[j];
            auto parentSet = state_->parentSetList[infection->id()];
            auto latentParent = state_->latentParents[j];

            sourceTransmissionProcessList.push_back(std::make_shared<SourceTransmissionImpl>(
                    coiProb,
                    state_->alleleFrequencies,
                    latentParent->loci(),
                    latentParent->latentGenotype()));

            transmissionProcessList.push_back(std::make_shared<TransmissionProcess>(
                    nodeTransmissionProcess,
                    sourceTransmissionProcessList.back(),
                    infection,
                    parentSet,
                    latentParent));
            likelihood.addTarget(transmissionProcessList.back());
        }
        this->setDirty();
    }

    Model::Model(State& state, double temperature) : Model(std::make_shared<State>(state), temperature) {}

    std::string Model::identifier() {
        return "ModelEight";
    }

    double Model::getTemperature() const {
        return temperature;
    }

    void Model::setTemperature(double t) {
        this->temperature = t;
        this->setDirty();
    }

    double Model::getPrior() {
        return prior.value();
    }

    double Model::getLikelihood() {
        return likelihood.value();
    }

    Likelihood Model::value() {
        if (isDirty()) {
            value_ = temperature * likelihood.value() + prior.value();
            this->setClean();
        }
        return value_;
    }
};// namespace transmission_nets::impl::ModelEight
