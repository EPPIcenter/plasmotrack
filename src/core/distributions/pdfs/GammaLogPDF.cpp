//
// Created by Maxwell Murphy on 7/27/20.
//


#include "GammaLogPDF.h"

#include <utility>


namespace transmission_nets::core::distributions {


    GammaLogPDF::GammaLogPDF(p_Parameterdouble target, p_Parameterdouble shape, p_Parameterdouble scale) : target_(std::move(target)), shape_(std::move(shape)), scale_(std::move(scale)) {
        target_->registerCacheableCheckpointTarget(this);
        target_->add_post_change_listener([=, this]() { this->setDirty(); });

        shape_->registerCacheableCheckpointTarget(this);
        shape_->add_post_change_listener([=, this]() { this->setDirty(); });

        scale_->registerCacheableCheckpointTarget(this);
        scale_->add_post_change_listener([=, this]() { this->setDirty(); });

        this->setDirty();
    }

    computation::Likelihood GammaLogPDF::value() {
        if (isDirty()) {
            logDenominator_ = std::lgamma(shape_->value()) + (shape_->value() * std::log(scale_->value()));
            value_          = (shape_->value() - 1) * std::log(target_->value()) +
                     (-target_->value() / scale_->value()) -
                     logDenominator_;
            this->setClean();
        }
        return value_;
    }

    std::string GammaLogPDF::identifier() {
        return {"GammaLogPDF"};
    }

}// namespace transmission_nets::core::distributions