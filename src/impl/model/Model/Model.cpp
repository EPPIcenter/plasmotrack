//
// Created by Maxwell Murphy on 11/15/23
//

#include "Model.h"
#include "ErrorRatePrior.h"

#include "core/abstract/observables/Dependency.h"
#include "core/distributions/pdfs/DiscretePDF.h"
#include "core/distributions/pdfs/GammaLogPDF.h"
#include "core/distributions/pdfs/BetaLogPDF.h"
#include "core/distributions/pdfs/FunctionLogPDF.h"

#include <cmath>
#include <limits>

#include <utility>

namespace transmission_nets::impl::Model {

    Model::Model(std::shared_ptr<State> state, double temperature) : state_(std::move(state)), temperature(temperature) {
        core::abstract::dependsOn(this, likelihood, prior);

        transmissionMemo        = std::make_shared<TransmissionMemoImpl>(state_->meanStrainsTransmitted, state_->alleleFrequencies, state_->transmissionLeak);
        if (state_->leakPrior) {
            // The memo follows the sampled leak.
            transmissionMemo->followLeak(state_->leakSize);
            prior.addTarget(std::make_shared<core::distributions::FunctionLogPDF>(
                    "transmission_leak_prior", [leak = state_->leakSize, p = *state_->leakPrior]() { return p.logDensity(leak->value()); }, state_->leakSize));
        }
        parentSetSizeLikelihood = std::make_shared<ParentSetSizeLikelihoodImpl>(state_->parentSetSizeProb);
        coiProb                 = std::make_shared<COIProbabilityImpl>(state_->meanCOI);

        // Register Priors
//        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->lossProb, state_->lossProbPriorAlpha, state_->lossProbPriorBeta));
        prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanCOI, state_->meanCOIPriorShape, state_->meanCOIPriorScale));
        prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->meanStrainsTransmitted, state_->meanStrainsTransmittedPriorShape, state_->meanStrainsTransmittedPriorScale));
        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->parentSetSizeProb, state_->parentSetSizePriorAlpha, state_->parentSetSizePriorBeta));
//        prior.addTarget(std::make_shared<core::distributions::BetaLogPDF>(state_->geometricGenerationProb, state_->geometricGenerationProbPriorAlpha, state_->geometricGenerationProbPriorBeta));
        if (state_->errorModel == ErrorModel::Shared) {
            // Every Infection holds the same two rates: one prior each.
            prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->expectedFalsePositives.front(), state_->obsFPRPriorShape, state_->obsFPRPriorScale));
            prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(state_->expectedFalseNegatives.front(), state_->obsFNRPriorShape, state_->obsFNRPriorScale));
        } else if (state_->errorModel == ErrorModel::Hierarchical) {
            // rate_i ~ Gamma(shape, mean / shape) truncated to the rates' range (ErrorRatePrior.h); mean ~ the wide Gamma
            // prior of its kind of rate; shape log-uniform on [HIERARCHICAL_SHAPE_LOWER, HIERARCHICAL_SHAPE_UPPER].
            const auto addRates = [&](const std::vector<State::p_Parameterdouble>& rates, const State::p_Parameterdouble& mean, const State::p_Parameterdouble& shape, double meanScale, const std::string& name) {
                for (const auto& rate : rates) {
                    prior.addTarget(std::make_shared<core::distributions::FunctionLogPDF>(
                            name, [=]() { return truncatedGammaLogDensity(rate->value(), mean->value(), shape->value()); }, rate, mean, shape));
                }
                prior.addTarget(std::make_shared<core::distributions::FunctionLogPDF>(
                        name + "_mean", [=]() { return core::distributions::gammaLogDensity(mean->value(), WIDE_RATE_PRIOR_SHAPE, meanScale); }, mean));
                prior.addTarget(std::make_shared<core::distributions::FunctionLogPDF>(
                        name + "_shape", [=]() {
                            const double s = shape->value();
                            return s >= HIERARCHICAL_SHAPE_LOWER and s <= HIERARCHICAL_SHAPE_UPPER ? -std::log(s) : -std::numeric_limits<double>::infinity();
                        },
                        shape));
            };
            addRates(state_->expectedFalsePositives, state_->efpMean, state_->efpShape, WIDE_FALSE_POSITIVE_PRIOR_SCALE, "efp");
            addRates(state_->expectedFalseNegatives, state_->efnMean, state_->efnShape, WIDE_FALSE_NEGATIVE_PRIOR_SCALE, "efn");
        } else {
            for (auto& obs : state_->expectedFalsePositives) {
                prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFPRPriorShape, state_->obsFPRPriorScale));
            }
            for (auto& obs : state_->expectedFalseNegatives) {
                prior.addTarget(std::make_shared<core::distributions::GammaLogPDF>(obs, state_->obsFNRPriorShape, state_->obsFNRPriorScale));
            }
        }

        int i = 0;
        for (auto& infection : state_->infections) {
            // infection duration likelihood

            if (!state_->nullModel) {
                if (infection->isSymptomatic()) {
                    prior.addTarget(std::make_shared<core::distributions::DiscretePDF<double>>(infection->infectionDuration(), state_->symptomaticInfectionDurationDist, "symptomatic_infection_duration"));
                } else {
                    prior.addTarget(std::make_shared<core::distributions::DiscretePDF<double>>(infection->infectionDuration(), state_->asymptomaticInfectionDurationDist, "asymptomatic_infection_duration"));
                }
            }

            for (auto& [locus, obsGenotype] : infection->observedGenotype()) {
                observationProcessLikelihoodList.push_back(std::make_shared<ObservationProcessImpl>(
                        obsGenotype,
                        infection->latentGenotype(locus),
                        state_->expectedFalsePositives[i],
                        state_->expectedFalseNegatives[i],
                        state_->nullModel
                        ));
                likelihood.addTarget(observationProcessLikelihoodList.back());
            }
            i++;
        }

            for (unsigned int j = 0; j < state_->infections.size(); ++j) {
                auto infection = state_->infections[j];
                // Names the infection's latent parent in its parent-set distribution; its genotypes are summed out.
                auto latentParent = state_->latentParents[j];

                transmissionProcessList.push_back(std::make_shared<TransmissionProcess>(
                        transmissionMemo,
                        coiProb,
                        parentSetSizeLikelihood,
                        infection,
                        state_->infectionEventOrdering,
                        latentParent,
                        state_->nullModel
                        ));
                likelihood.addTarget(transmissionProcessList.back());
            }

        this->setDirty();
    }

    Model::Model(State& state, double temperature) : Model(std::make_shared<State>(state), temperature) {}

    std::string Model::identifier() {
        return "Model";
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
};// namespace transmission_nets::impl::Model
