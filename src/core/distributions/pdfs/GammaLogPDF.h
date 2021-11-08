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

    class GammaLogPDF : public computation::PartialLikelihood  {
    public:
        using p_ParameterDouble = std::shared_ptr<core::parameters::Parameter<double>>;
        GammaLogPDF(p_ParameterDouble target, p_ParameterDouble shape, p_ParameterDouble scale);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_ParameterDouble target_;
        p_ParameterDouble shape_;
        p_ParameterDouble scale_;
        double logDenominator_ = 0;
    };

}



#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
