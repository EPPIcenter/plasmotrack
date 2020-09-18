//
// Created by Maxwell Murphy on 6/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H
#define TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H

#include <vector>
#include <map>

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/parameters/TransmissionNetwork.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Parameter.h"


namespace transmission_nets::impl {

    struct ModelThreeState {
        static constexpr int MAX_ALLELES = 32;
        static constexpr int MAX_PARENT_SET = 1;
        using LocusImpl = core::containers::Locus;

        using GeneticsImpl = core::datatypes::AllelesBitSet<MAX_ALLELES>;
        using InfectionEvent = core::containers::Infection<GeneticsImpl, LocusImpl>;

        using AlleleFrequencyImpl = core::datatypes::Simplex;
        using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;


        std::map<std::string, LocusImpl *> loci{};
        std::vector<InfectionEvent *> infections{};

        AlleleFrequencyContainerImpl alleleFrequencies;

        // Network Structure
        core::parameters::TransmissionNetwork<InfectionEvent> transmissionNetwork;

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




#endif //TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H
