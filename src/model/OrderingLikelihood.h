//
// Created by Maxwell Murphy on 2/25/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERINGLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_ORDERINGLIKELIHOOD_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/Ordering.h"

namespace transmission_nets::model::ordering {

    template<typename TimeDiffCDF, typename ElementType>
    class OrderingLikelihood : public core::computation::PartialLikelihood {

    public:
        OrderingLikelihood(core::parameters::Ordering<ElementType>& ordering, TimeDiffCDF& tdiff);

        core::computation::Likelihood value() override;
        std::string identifier() override;

    private:
        core::parameters::Ordering<ElementType>& ordering_;
        TimeDiffCDF& tdiff_;
    };


    template<typename TimeDiffCDF, typename ElementType>
    OrderingLikelihood<TimeDiffCDF, ElementType>::OrderingLikelihood(core::parameters::Ordering<ElementType>& ordering, TimeDiffCDF& tdiff) : ordering_(Ordering), tdiff_(tdiff) {
    }

    template<typename TimeDiffCDF, typename ElementType>
    core::computation::Likelihood OrderingLikelihood<TimeDiffCDF, ElementType>::value() {
        //Todo
        return PartialLikelihood::value();
    }

    template<typename TimeDiffCDF, typename ElementType>
    std::string OrderingLikelihood<TimeDiffCDF, ElementType>::identifier() {
        return std::string("OrderingLikelihood");
    }

}// namespace transmission_nets::model::ordering

#endif//TRANSMISSION_NETWORKS_APP_ORDERINGLIKELIHOOD_H
