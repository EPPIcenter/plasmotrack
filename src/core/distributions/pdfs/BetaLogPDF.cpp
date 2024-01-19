//
// Created by Maxwell Murphy on 7/27/20.
//


#include "BetaLogPDF.h"

#include <utility>

namespace transmission_nets::core::distributions {

    BetaLogPDF::BetaLogPDF(p_Parameterdouble target, p_Parameterdouble alpha, p_Parameterdouble beta) : target_(std::move(target)), alpha_(std::move(alpha)), beta_(std::move(beta)) {
        target_->registerCacheableCheckpointTarget(this);
        target_->add_post_change_listener([=, this]() { this->setDirty(); });

        alpha_->registerCacheableCheckpointTarget(this);
        alpha_->add_post_change_listener([=, this]() { this->setDirty(); });

        beta_->registerCacheableCheckpointTarget(this);
        beta_->add_post_change_listener([=, this]() { this->setDirty(); });

        this->setDirty();
    }


    computation::Likelihood BetaLogPDF::value() {
        if (isDirty()) {
            logDenominator_ = std::lgamma(alpha_->value()) + std::lgamma(beta_->value()) - std::lgamma(alpha_->value() + beta_->value());
            value_          = (alpha_->value() - 1) * log(target_->value()) +
                     (beta_->value() - 1) * log(1 - target_->value()) -
                     logDenominator_;
            setClean();
        }

        return value_;
    }

    std::string BetaLogPDF::identifier() {
        return {"BetaLogPDF"};
    }

}// namespace transmission_nets::core::distributions
