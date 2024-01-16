//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_BETALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

#include <cmath>
#include <memory>

namespace transmission_nets::core::distributions {

    class BetaLogPDF : public computation::PartialLikelihood {
    public:
        using p_Parameterfloat = std::shared_ptr<core::parameters::Parameter<float>>;
        BetaLogPDF(p_Parameterfloat target, p_Parameterfloat alpha, p_Parameterfloat beta);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_Parameterfloat target_;
        p_Parameterfloat alpha_;
        p_Parameterfloat beta_;
        float logDenominator_ = 0;
    };

}// namespace transmission_nets::core::distributions

#endif//TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
