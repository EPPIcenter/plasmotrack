//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_BETALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

class BetaLogPDF : public PartialLikelihood {
public:
    BetaLogPDF(Parameter<double> &target, double alpha, double beta);
    Likelihood value() override;

private:
    Parameter<double> &target_;
    double alpha_;
    double beta_;
    double logDenominator_;
};



#endif//TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
