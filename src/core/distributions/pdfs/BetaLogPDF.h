//
// Created by Maxwell Murphy on 7/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
#define TRANSMISSION_NETWORKS_APP_BETALOGPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Parameter.h"

namespace transmission_nets::core::distributions {

    class BetaLogPDF : public computation::PartialLikelihood {
    public:
        BetaLogPDF(parameters::Parameter<double> &target, double alpha, double beta);
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        parameters::Parameter<double> &target_;
        double alpha_;
        double beta_;
        double logDenominator_;
    };

}

#endif//TRANSMISSION_NETWORKS_APP_BETALOGPDF_H
