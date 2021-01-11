//
// Created by Maxwell Murphy on 6/2/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NETWORKBASEDTRANSMISSIONPROCESS_H
#define TRANSMISSION_NETWORKS_APP_NETWORKBASEDTRANSMISSIONPROCESS_H

#include "core/computation/PartialLikelihood.h"
#include "core/parameters/TransmissionNetwork.h"

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
    class NetworkBasedTransmissionProcess : public core::computation::PartialLikelihood {

    public:
        NetworkBasedTransmissionProcess(NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp,
                                        InfectionEventImpl &child, core::parameters::Parameter<core::containers::ParentSet<InfectionEventImpl>> &parentSet);

        Likelihood value() override;
        std::string identifier() override;

    private:
        NodeTransmissionProcessImpl &ntp_;
        SourceTransmissionProcessImpl &stp_;
        InfectionEventImpl &child_;
        core::parameters::Parameter<core::containers::ParentSet<InfectionEventImpl>> &parentSet_;


    };

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
    std::string NetworkBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::identifier() {
        return "NetworkBasedTransmissionProcess";
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
    Likelihood
    NetworkBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::value() {
        if (isDirty()) {
            if(parentSet_.value().size() == 0) {
                this->value_ = stp_.value();
            } else {
                this->value_ = ntp_.calculateLogLikelihood(child_, parentSet_.value());
            }
            setClean();
        }

        return value_;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl>
    NetworkBasedTransmissionProcess<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl>::NetworkBasedTransmissionProcess(
            NodeTransmissionProcessImpl &ntp, SourceTransmissionProcessImpl &stp, InfectionEventImpl &child,
            core::parameters::Parameter<core::containers::ParentSet<InfectionEventImpl>> &parentSet):ntp_(ntp), stp_(stp), child_(child), parentSet_(parentSet) {
        ntp_.add_set_dirty_listener([=, this]() { setDirty(); });
        ntp_.registerCacheableCheckpointTarget(this);

        stp_.add_set_dirty_listener([=, this]() { setDirty(); });
        stp_.registerCacheableCheckpointTarget(this);

        child_.add_post_change_listener([=, this]() { setDirty(); });
        child_.registerCacheableCheckpointTarget(this);

        parentSet_.add_post_change_listener([=, this]() { setDirty(); });
        parentSet_.registerCacheableCheckpointTarget(this);
    }

}



#endif //TRANSMISSION_NETWORKS_APP_NETWORKBASEDTRANSMISSIONPROCESS_H
