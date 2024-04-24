//
// Created by Maxwell Murphy on 11/15/23
//

#include "Model.h"

#include <utility>

namespace transmission_nets::impl::ModelNine {

    Model::Model(std::shared_ptr<State> state, double temperature) : state_(std::move(state)), temperature(temperature) {
        likelihood.add_set_dirty_listener([=, this]() {
            this->setDirty();
        });
        likelihood.registerCacheableCheckpointTarget(this);

        prior.add_set_dirty_listener([=, this]() {
            this->setDirty();
        });
        prior.registerCacheableCheckpointTarget(this);

//        intp                    = std::make_shared<InterTransmissionProbImpl>(state_->geometricGenerationProb);
//        nodeTransmissionProcess = std::make_shared<NodeTransmissionImpl>(state_->lossProb, intp);
        nodeTransmissionProcess = std::make_shared<NodeTransmissionImpl>(state_->meanStrainsTransmitted);
        coiProb                 = std::make_shared<COIProbabilityImpl>(state_->meanCOI);

        // Register Priors
//        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->lossProb, state_->lossProbPriorAlpha, state_->lossProbPriorBeta));
        prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanCOI, state_->meanCOIPriorShape, state_->meanCOIPriorScale));
        prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanStrainsTransmitted, state_->meanStrainsTransmittedPriorShape, state_->meanStrainsTransmittedPriorScale));
//        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->geometricGenerationProb, state_->geometricGenerationProbPriorAlpha, state_->geometricGenerationProbPriorBeta));
        for (auto& obs : state_->expectedFalsePositives) {
            prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFPRPriorShape, state_->obsFPRPriorScale));
        }
        for (auto& obs : state_->expectedFalseNegatives) {
            prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFNRPriorShape, state_->obsFNRPriorScale));
        }

        int i = 0;
        for (auto& infection : state_->infections) {
            // infection duration likelihood

            if (infection->isSymptomatic()) {
                prior.addTarget(std::make_shared<core::distributions::DiscretePDF<double>>(infection->infectionDuration(), state_->symptomaticInfectionDurationDist, "symptomatic_infection_duration"));
            } else {
                prior.addTarget(std::make_shared<core::distributions::DiscretePDF<double>>(infection->infectionDuration(), state_->asymptomaticInfectionDurationDist, "asymptomatic_infection_duration"));
            }


                for (auto& [locus, obsGenotype] : infection->observedGenotype()) {
                    observationProcessLikelihoodList.push_back(std::make_shared<ObservationProcessImpl>(
                            obsGenotype,
                            infection->latentGenotype(locus),
                            state_->expectedFalsePositives[i],
                            state_->expectedFalseNegatives[i],
                            state_->null_model_));
                    likelihood.addTarget(observationProcessLikelihoodList.back());
                }

            state_->parentSetList[infection->id()] = std::make_shared<ParentSetImpl>(state_->infectionEventOrdering, infection, state_->allowedRelationships->allowedParents(infection));
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
                        latentParent->latentGenotype(),
                        state_->null_model_));

                transmissionProcessList.push_back(std::make_shared<TransmissionProcess>(
                        nodeTransmissionProcess,
                        sourceTransmissionProcessList.back(),
                        infection,
                        parentSet,
                        latentParent,
                        state_->null_model_));
                likelihood.addTarget(transmissionProcessList.back());
            }

        this->setDirty();
    }

    Model::Model(State& state, double temperature) : Model(std::make_shared<State>(state), temperature) {}

    std::string Model::identifier() {
        return "ModelNine";
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
            std::unique_lock lock(valueMutex);
            value_ = temperature * likelihood.value() + prior.value();
            lock.unlock();
            this->setClean();
        }
        return value_;
    }

    Likelihood Model::valueThreadSafe() {
        std::shared_lock lock(valueMutex);
        auto val = value_;
        lock.unlock();
        return val;
    }
};// namespace transmission_nets::impl::ModelNine
