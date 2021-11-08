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
        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
        BetaLogPDF(p_ParameterDouble target, p_ParameterDouble alpha, p_ParameterDouble beta);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_ParameterDouble target_;
        p_ParameterDouble alpha_;
        p_ParameterDouble beta_;
        double logDenominator_ = 0;
    };

}

#endif//TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
