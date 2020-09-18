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
        GammaLogPDF(parameters::Parameter<double> &target, double shape, double scale);
        computation::Likelihood value() override;

    private:
        parameters::Parameter<double> &target_;
        double shape_;
        double scale_;
        double logDenominator_;
    };

}



#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
