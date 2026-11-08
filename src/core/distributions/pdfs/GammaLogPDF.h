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

    // log of the Gamma(shape, scale) density at x.
    inline double gammaLogDensity(const double x, const double shape, const double scale) {
        const double logDenominator = std::lgamma(shape) + (shape * std::log(scale));
        return (shape - 1) * std::log(x) + (-x / scale) - logDenominator;
    }

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
    };

}// namespace transmission_nets::core::distributions


#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
