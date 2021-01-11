//
// Created by Maxwell Murphy on 3/10/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MUTABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_MUTABLEPASSTHROUGH_H

#include <tuple>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {

    template<typename T>
    class MutablePassthrough : public crtp<T, MutablePassthrough> {

        using ChangeCallback = std::function<void()>;
        CREATE_EVENT(pre_change, ChangeCallback)
        CREATE_EVENT(post_change, ChangeCallback)

    public:

        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0& target);

        void saveState() noexcept;

        void restoreState() noexcept;

        void acceptState() noexcept;

    };

    template<typename T>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> MutablePassthrough<T>::registerCheckpointTarget(T0 &target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=, this]() { target.saveState(); });
        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=, this]() { target.acceptState(); });
        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=, this]() { target.restoreState(); });
        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename T>
    void MutablePassthrough<T>::saveState() noexcept {
        this->underlying().notify_save_state();
    }

    template<typename T>
    void MutablePassthrough<T>::restoreState() noexcept {
        this->underlying().notify_restore_state();
    }

    template<typename T>
    void MutablePassthrough<T>::acceptState() noexcept {
        this->underlying().notify_accept_state();
    }

}


#endif //TRANSMISSION_NETWORKS_APP_MUTABLEPASSTHROUGH_H
