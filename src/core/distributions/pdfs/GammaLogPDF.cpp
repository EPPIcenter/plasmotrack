//
// Created by Maxwell Murphy on 7/27/20.
//

#include <cmath>

#include "GammaLogPDF.h"


namespace transmission_nets::core::distributions {

    GammaLogPDF::GammaLogPDF(parameters::Parameter<double> &target, double shape, double scale) :  target_(target), shape_(shape), scale_(scale) {
        target_.registerCacheableCheckpointTarget(this);
        target_.add_post_change_listener([=, this]() { this->setDirty(); });

        logDenominator_ = lgamma(shape_) - (shape_ * log(scale_));
        this->setDirty();
    }

    computation::Likelihood GammaLogPDF::value() {
        if(isDirty()) {
            value_ = (shape_ - 1) * log(target_.value()) +
                     (-target_.value() / scale_) -
                     logDenominator_;
            this->setClean();
        }
        return value_;
    }

}