//
// Created by Maxwell Murphy on 7/27/20.
//


#include "GammaLogPDF.h"
#include "core/abstract/observables/Dependency.h"

#include <utility>


namespace transmission_nets::core::distributions {


    GammaLogPDF::GammaLogPDF(p_Parameterdouble target, p_Parameterdouble shape, p_Parameterdouble scale) : target_(std::move(target)), shape_(std::move(shape)), scale_(std::move(scale)) {
        abstract::dependsOn(this, target_, shape_, scale_);
        this->setDirty();
    }

    computation::Likelihood GammaLogPDF::value() {
        if (isDirty()) {
            value_ = gammaLogDensity(target_->value(), shape_->value(), scale_->value());
            this->setClean();
        }
        return value_;
    }

    std::string GammaLogPDF::identifier() {
        return {"GammaLogPDF"};
    }

}// namespace transmission_nets::core::distributions