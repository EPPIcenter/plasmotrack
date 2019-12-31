//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H


#include "LikelihoodFunction.h"

typedef float LogLik;

template <class LikelihoodFunctor, class LikelihoodFunctorParameters>
class PartialLikelihood {
public:
    PartialLikelihood(LikelihoodFunctor likelihoodFunctor, LikelihoodFunctorParameters likelihoodFunctorParameters);
    LogLik get_value();
    void rollback();
    void invalidate_cache();

private:
    bool cache_valid { false };
    LogLik cached_value { 0 };
    LogLik previous_value { 0 };
    LikelihoodFunctor likelihood_functor;
    LikelihoodFunctorParameters likelihood_functor_parameters;
};

template<class LikelihoodFunctor, class LikelihoodFunctorParameters>
LogLik PartialLikelihood<LikelihoodFunctor, LikelihoodFunctorParameters>::get_value() {
    if (!cache_valid) {
        previous_value = cached_value;
        cached_value = likelihood_functor(likelihood_functor_parameters);
        cache_valid = true;
    }
    return cached_value;
}

template<class LikelihoodFunctor, class LikelihoodFunctorParameters>
void PartialLikelihood<LikelihoodFunctor, LikelihoodFunctorParameters>::invalidate_cache() {
    cache_valid = false;
}

template<class LikelihoodFunctor, class LikelihoodFunctorParameters>
void PartialLikelihood<LikelihoodFunctor, LikelihoodFunctorParameters>::rollback() {
    cached_value = previous_value;
}

template<class LikelihoodFunctor, class LikelihoodFunctorParameters>
PartialLikelihood<LikelihoodFunctor, LikelihoodFunctorParameters>::PartialLikelihood(
        LikelihoodFunctor likelihoodFunctor, LikelihoodFunctorParameters likelihoodFunctorParameters)
        :likelihood_functor(likelihoodFunctor), likelihood_functor_parameters(likelihoodFunctorParameters) {
            likelihoodFunctorParameters.register_likelihood(this);
        }


#endif //TRANSMISSION_NETWORKS_APP_PARTIALLIKELIHOOD_H