//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

class GammaLogPDF : public PartialLikelihood  {
public:
    GammaLogPDF(Parameter<double> &target, double shape, double scale);
    Likelihood value() override;

private:
    Parameter<double> &target_;
    double shape_;
    double scale_;
    double logDenominator_;
};


#endif//TRANSMISSION_NETWORKS_APP_GAMMALOGPDF_H
