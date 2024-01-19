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
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        GammaLogPDF(p_Parameterdouble target, p_Parameterdouble shape, p_Parameterdouble scale);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_Parameterdouble target_;
        p_Parameterdouble shape_;
        p_Parameterdouble scale_;
        double logDenominator_ = 0;
    };

}// namespace transmission_nets::core::distributions


#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
