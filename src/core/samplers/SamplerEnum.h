//
// Created by Maxwell Murphy on 4/10/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H
#define TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H

namespace transmission_nets::core::samplers {
    enum SAMPLER_STATE_ID {
        ContinuousRandomWalkID,
        ConstrainedContinuousRandomWalkID,
        SALTID,
        RandomAlleleBitSet3ID,
        RandomAlleleBitSet4ID,
        ParentRelinkID,
        CloneOrderSwapID,
        InformedOrderJumpID,
        InformedOrderJumpTryID,
        ParentRelinkTryID,
        AlleleSwapID,
        StrainCountID
    };
}

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERENUM_H
