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
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        BetaLogPDF(p_Parameterdouble target, p_Parameterdouble alpha, p_Parameterdouble beta);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_Parameterdouble target_;
        p_Parameterdouble alpha_;
        p_Parameterdouble beta_;
        double logDenominator_ = 0;
    };

}// namespace transmission_nets::core::distributions

#endif//TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
