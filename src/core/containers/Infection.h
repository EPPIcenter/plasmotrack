//
// Created by Maxwell Murphy on 1/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFECTION_H
#define TRANSMISSION_NETWORKS_APP_INFECTION_H

#include <boost/container/flat_map.hpp>

#include "core/containers/Locus.h"
#include "core/datatypes/Data.h"
#include "core/parameters/Parameter.h"

template<template <typename> typename Wrapper, typename GeneticImpl, typename LocusImpl>
using GenotypeMap = boost::container::flat_map<LocusImpl*, Wrapper<GeneticImpl>>; // Maps genetic readout to loci

template<typename GeneticImpl, typename LocusImpl = Locus>
using LocusAssignment = std::pair<LocusImpl*, GeneticImpl>;

template <typename GeneticImpl, typename LocusImpl = Locus>
class Infection: public Observable<Infection<GeneticImpl, LocusImpl>> {
public:

    template<typename LocusDataIter>
    Infection(LocusDataIter obs, LocusDataIter latent) {
        for(auto const& [locus, genetics] : obs) {
            observedGenotype_.emplace(locus, genetics);
        }

        for(auto const& [locus, genetics] : latent) {
            latentGenotype_.emplace(locus, genetics);
        }
    };

    GenotypeMap<Data, GeneticImpl, LocusImpl>& observedGenotype() {
        return observedGenotype_;
    };

    Data<GeneticImpl>& observedGenotype(LocusImpl* locus) {
        return observedGenotype_.at(locus);
    };

    GenotypeMap<Parameter, GeneticImpl, LocusImpl>& latentGenotype() {
        return latentGenotype_;
    };

    Parameter<GeneticImpl>& latentGenotype(LocusImpl* locus) {
        return latentGenotype_.at(locus);
    };

private:
    GenotypeMap<Data, GeneticImpl, LocusImpl> observedGenotype_{};
    GenotypeMap<Parameter, GeneticImpl, LocusImpl> latentGenotype_{};
};

#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
