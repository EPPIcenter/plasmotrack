//
// Created by Maxwell Murphy on 6/5/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELTWOSTATE_H
#define TRANSMISSION_NETWORKS_APP_MODELTWOSTATE_H

#include <vector>
#include <map>

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Parameter.h"
#include "core/parameters/Ordering.h"


namespace transmission_nets::impl {

    struct ModelTwoState {
        static constexpr int MAX_ALLELES = 32;
        using LocusImpl = core::containers::Locus;
        using GeneticsImpl = core::datatypes::AllelesBitSet<MAX_ALLELES>;
        using InfectionEvent = core::containers::Infection<GeneticsImpl, LocusImpl>;
        using AlleleFrequencyImpl = core::datatypes::Simplex;
        using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;


        std::map<std::string, LocusImpl *> loci{};
        std::vector<InfectionEvent *> infections{};

        AlleleFrequencyContainerImpl alleleFrequencies;

        // Network Structure
        core::parameters::Ordering<InfectionEvent> infectionEventOrdering;

        // Observation Process
        core::parameters::Parameter<double> observationFalsePositiveRate;
        core::parameters::Parameter<double> observationFalseNegativeRate;

        // Node Transmission Process
        core::parameters::Parameter<double> geometricGenerationProb;
        core::parameters::Parameter<double> lossProb;
        core::parameters::Parameter<double> mutationProb;

        // Source Transmission Process
        core::parameters::Parameter<double> meanCOI;
    };

}




#endif //TRANSMISSION_NETWORKS_APP_MODELTWOSTATE_H
