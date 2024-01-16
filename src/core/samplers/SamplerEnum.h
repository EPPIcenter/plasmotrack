//
// Created by Maxwell Murphy on 4/10/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H
#define TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H

namespace transmission_nets::core::samplers {
    enum SAMPLER_STATE_ID {
        ContinuousRandomWalkID,
        ConstrainedContinuousRandomWalkID,
        DiscreteRandomWalkID,
        ConstrainedDiscreteRandomWalkID,
        SALTID,
        SimplexID,
        RandomAlleleBitSetID,
        RandomAlleleBitSet2ID,
        RandomAlleleBitSet3ID,
        RandomAlleleBitSet4ID,
        SequentialAllelesBitSetID,
        ZanellaAllelesBitSetID,
        ZanellaAllelesBitSetFlipID,
        JointGeneticsTimeID,
        ZanellaJointGeneticsGraphID,
        OrderID,
        RandomAddEdgeID,
        RandomRemoveEdgeID,
        RandomReverseEdgeID,
        RandomSwapEdgeID,
        ZanellaNeighborOrderID,
        ZanellaOrderID
    };
}

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H
