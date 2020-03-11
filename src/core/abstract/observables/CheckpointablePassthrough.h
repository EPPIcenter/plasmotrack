//
// Created by Maxwell Murphy on 3/10/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H


#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

template<typename T>
class CheckpointablePassthrough : public crtp<T, CheckpointablePassthrough> {
    using CallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(save_state, CallbackType);
    CRTP_CREATE_EVENT(accept_state, CallbackType);
    CRTP_CREATE_EVENT(restore_state, CallbackType);

public:

    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0& target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([&]() { target.saveState(); });
        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([&]() { target.acceptState(); });
        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([&]() { target.restoreState(); });
        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    void saveState() noexcept {
        this->underlying().notify_save_state();
    }

    void restoreState() noexcept {
        this->underlying().notify_restore_state();
    }

    void acceptState() noexcept {
        this->underlying().notify_accept_state();
    }

//    bool constexpr isSaved() noexcept {
//        return true;
//    }

};

#endif //TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
