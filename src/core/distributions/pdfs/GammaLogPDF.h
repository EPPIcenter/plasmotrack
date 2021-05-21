//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

namespace transmission_nets::core::distributions {

    class GammaLogPDF : public computation::PartialLikelihood  {
    public:
        GammaLogPDF(parameters::Parameter<double> &target, parameters::Parameter<double> &shape, parameters::Parameter<double> &scale);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        parameters::Parameter<double> &target_;
        parameters::Parameter<double> &shape_;
        parameters::Parameter<double> &scale_;
        double logDenominator_;
    };

}



#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
