//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/containers/Infection.h"

#include "core/computation/Computation.h"

template<typename COIProbabilityImpl, typename AlleleFrequencyContainer, typename GeneticsImpl>
class MultinomialSourceTransmissionProcess : public Computation<double>,
                                             public Observable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GeneticsImpl>>,
                                             public Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GeneticsImpl>>,
                                             public Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GeneticsImpl>, double> {
    using CallbackType = std::function<void()>;
    CREATE_EVENT(save_state, CallbackType)
    CREATE_EVENT(accept_state, CallbackType)
    CREATE_EVENT(restore_state, CallbackType)

public:

    MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb,
                                         AlleleFrequencyContainer &alleleFrequenciesContainer,
                                         Infection<GeneticsImpl> &founder)
            : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer), founder_(founder) {

        coiProb_.registerCacheableCheckpointTarget(*this);
        coiProb_.registerDirtyTarget(*this);

        alleleFrequenciesContainer_.registerCacheableCheckpointTarget(*this);
        alleleFrequenciesContainer_.add_post_change_listener([&]() { this->setDirty(); });

        founder_.registerCacheableCheckpointTarget(*this);
        founder_.add_post_change_listener([&]() { this->setDirty(); });
    }

    double value() override {
        if (this->isDirty()) {
            this->value_ = calculateLikelihood();
            this->setClean();
        }
        return this->value_;
    };


private:
    friend class Cacheable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GeneticsImpl>>;
    friend class Checkpointable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, GeneticsImpl>, double>;

    COIProbabilityImpl &coiProb_;
    AlleleFrequencyContainer &alleleFrequenciesContainer_;
    Infection<GeneticsImpl> &founder_;

    double calculateLogLikelihood() {
        double llik = 0.0;
        const auto founderLatentGenotype = founder_.latentGenotype();
        for (auto const& [locus, founderGenotypeAtLocus] : founderLatentGenotype) {
            double probNotObserved = 0.0;
            const double totalAlleles = locus->totalAlleles();
            const auto alleleFrequenciesParameter = alleleFrequenciesContainer_.alleleFrequencies(locus);
            // Calculate prob that allele was not drawn in one infection
            for (int j = 0; j < totalAlleles; ++j) {
                if(!founderGenotypeAtLocus.value().allele(j)) {
                    probNotObserved += alleleFrequenciesParameter.value().frequencies(j);
                }
            }

            // Calculate likelihood that allele was not drawn in k infections times prob of k infections
            double locusLik = 0.0;
            unsigned int minCOI = founderGenotypeAtLocus.value().totalPositiveCount();
            unsigned int totalInfections = coiProb_.value().size();
            for (auto k = minCOI; k < totalInfections; ++k) {
                locusLik += pow(probNotObserved, k) * coiProb_.value()(k);
            }
            llik += log(1 - locusLik);
        }
        return llik;
    }

    double calculateLikelihood() {
        return exp(calculateLogLikelihood());
    }
};

#endif //TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
