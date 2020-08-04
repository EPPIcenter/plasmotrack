//
// Created by Maxwell Murphy on 7/27/20.
//

#include <cmath>

#include "BetaLogPDF.h"

BetaLogPDF::BetaLogPDF(Parameter<double> &target, double alpha, double beta) : target_(target), alpha_(alpha), beta_(beta) {
    target_.registerCacheableCheckpointTarget(this);
    target_.add_post_change_listener([=]() { this->setDirty(); });

    logDenominator_ = lgamma(alpha_) + lgamma(beta_) - lgamma(alpha_ + beta_);

    this->setDirty();
}


Likelihood BetaLogPDF::value() {
    if (isDirty()) {
        value_ = (alpha_ - 1) * log(target_.value()) +
                 (beta_ - 1) * log(1 - target_.value()) -
                 logDenominator_;
        setClean();
    }
    return value_;
}
