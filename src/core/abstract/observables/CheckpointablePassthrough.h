//
// CheckpointablePassthrough: a CRTP mixin for a container of Checkpointable parts, which announces their saves,
// accepts and restores as its own (see passThrough in Dependency.h). It holds no saved states itself.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H


#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

#include <functional>

namespace transmission_nets::core::abstract {

    template<typename T>
    class CheckpointablePassthrough : public crtp<T, CheckpointablePassthrough> {
        using AcceptCallbackType      = std::function<void()>;
        using SaveRestoreCallbackType = std::function<void(int saved_state_id)>;
        CRTP_CREATE_EVENT(save_state, SaveRestoreCallbackType)
        CRTP_CREATE_EVENT(accept_state, AcceptCallbackType)
        CRTP_CREATE_EVENT(restore_state, SaveRestoreCallbackType)
    };

}// namespace transmission_nets::core::abstract


#endif//TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
