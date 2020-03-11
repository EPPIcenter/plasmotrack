//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/CacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"

template<typename COIProbabilityImpl, typename AlleleFrequencyImpl>
class MultinomialSourceTransmissionProcess : public Observable<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyImpl>>,
                                             public CacheablePassthrough<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyImpl>>,
                                             public CheckpointablePassthrough<MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyImpl>> {
    using CallbackType = std::function<void()>;
    CREATE_EVENT(save_state, CallbackType);
    CREATE_EVENT(accept_state, CallbackType);
    CREATE_EVENT(restore_state, CallbackType);

public:

    MultinomialSourceTransmissionProcess(COIProbabilityImpl &coiProb,
                                         AlleleFrequencyContainer<AlleleFrequencyImpl> &alleleFrequenciesContainer)
            : coiProb_(coiProb), alleleFrequenciesContainer_(alleleFrequenciesContainer) {

        coiProb_.registerCheckpointTarget(*this);
        coiProb_.registerDirtyTarget(*this);
        alleleFrequenciesContainer_.registerCheckpointTarget(*this);
        alleleFrequenciesContainer_.registerDirtyTarget(*this);
    }

    template<typename GeneticsImpl>
    double calculateLogLikelihood(Infection<GeneticsImpl>& founder) {
        double llik = 0.0;
        const auto founderLatentGenotype = founder.latentGenotype();
        for (auto const& [locus, founderGenotypeAtLocus] : founderLatentGenotype) {
            double probNotObserved = 0.0;
            const double totalAlleles = locus->totalAlleles();
            const auto alleleFrequenciesParameter = alleleFrequenciesContainer_.alleleFrequencies(locus);
            // Calculate prob that allele was not drawn in one infection
            for (int j = 0; j < totalAlleles; ++j) {
                if(!founderGenotypeAtLocus.value().allele(j)) {
                    probNotObserved += alleleFrequenciesParameter.value().alleleFrequencies(j);
                }
            }

            // Calculate likelihood that allele was not drawn in k infections times prob of k infections
            double locusLik = 0.0;
            size_t totalInfections = coiProb_.value().size();
            for (size_t k = 1; k < totalInfections; ++k) {
                locusLik += pow(probNotObserved, k) * coiProb_.value()(k);
            }
            llik += log(locusLik);
        }
        return llik;
    }

    template<typename GeneticsImpl>
    double calculateLikelihood(Infection<GeneticsImpl>& founder) {
        return exp(calculateLogLikelihood(founder));
    }


private:
    COIProbabilityImpl &coiProb_;
    AlleleFrequencyContainer<AlleleFrequencyImpl> &alleleFrequenciesContainer_;
};

#endif //TRANSMISSION_NETWORKS_APP_MULTINOMIALSOURCETRANSMISSIONPROCESS_H
