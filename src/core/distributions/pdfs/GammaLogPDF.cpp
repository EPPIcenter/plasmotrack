//
// Created by Maxwell Murphy on 7/27/20.
//

#include <cmath>

#include "GammaLogPDF.h"


namespace transmission_nets::core::distributions {


    GammaLogPDF::GammaLogPDF(parameters::Parameter<double> &target, parameters::Parameter<double> &shape, parameters::Parameter<double> &scale) :  target_(target), shape_(shape), scale_(scale) {
        target_.registerCacheableCheckpointTarget(this);
        target_.add_post_change_listener([=, this]() { this->setDirty(); });

        shape_.registerCacheableCheckpointTarget(this);
        shape_.add_post_change_listener([=, this]() { this->setDirty(); });

        scale_.registerCacheableCheckpointTarget(this);
        scale_.add_post_change_listener([=, this]() { this->setDirty(); });

        this->setDirty();
    }

    computation::Likelihood GammaLogPDF::value() {
        if(isDirty()) {
            logDenominator_ = std::lgamma(shape_.value()) + (shape_.value() * std::log(scale_.value()));
            value_ = (shape_.value() - 1) * log(target_.value()) +
                     (-target_.value() / scale_.value()) -
                     logDenominator_;
            this->setClean();
        }
        return value_;
    }

    std::string GammaLogPDF::identifier() {
        return std::string("GammaLogPDF");
    }

}