//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

#include <cmath>
#include <memory>

namespace transmission_nets::core::distributions {

    class GammaLogPDF : public computation::PartialLikelihood {
    public:
        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;
        GammaLogPDF(p_Parameterfloat target, p_Parameterfloat shape, p_Parameterfloat scale);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_Parameterfloat target_;
        p_Parameterfloat shape_;
        p_Parameterfloat scale_;
        float logDenominator_ = 0;
    };

}// namespace transmission_nets::core::distributions


#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
