//
// Created by Maxwell Murphy on 3/10/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H


#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {

    template<typename T>
    class CheckpointablePassthrough : public crtp<T, CheckpointablePassthrough> {
        using CallbackType = std::function<void()>;
        using AcceptCallbackType = std::function<void()>;
        using SaveRestoreCallbackType = std::function<void(const std::string& saved_state_id)>;
        CRTP_CREATE_EVENT(save_state, SaveRestoreCallbackType)
        CRTP_CREATE_EVENT(accept_state, AcceptCallbackType)
        CRTP_CREATE_EVENT(restore_state, SaveRestoreCallbackType)

        struct StateCheckpoint {
            explicit StateCheckpoint(std::string savedStateId) : saved_state_id(std::move(savedStateId)) {}
            std::string saved_state_id;
        };


    public:

        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0 *target);

        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCacheableCheckpointTarget(T0 *target);

        void addPreSaveHook(const CallbackType& cb) {
            this->pre_save_hooks_.emplace_back(cb);
        }

        void addPostSaveHook(const CallbackType& cb) {
            this->post_save_hooks_.emplace_back(cb);
        }

        void addPreRestoreHook(const CallbackType& cb) {
            this->pre_restore_hooks_.emplace_back(cb);
        }

        void addPostRestoreHook(const CallbackType& cb) {
            this->post_restore_hooks_.emplace_back(cb);
        }

        void addPreAcceptHook(const CallbackType& cb) {
            this->pre_accept_hooks_.emplace_back(cb);
        }

        void addPostAcceptHook(const CallbackType& cb) {
            this->post_accept_hooks_.emplace_back(cb);
        }

        void saveState(std::string savedStateId) noexcept;

        void restoreState(std::string savedStateId) noexcept;

        void acceptState() noexcept;

    protected:
        std::vector<StateCheckpoint> saved_states_stack_{};
        std::vector<CallbackType> pre_save_hooks_{};
        std::vector<CallbackType> post_save_hooks_{};
        std::vector<CallbackType> pre_restore_hooks_{};
        std::vector<CallbackType> post_restore_hooks_{};
        std::vector<CallbackType> pre_accept_hooks_{};
        std::vector<CallbackType> post_accept_hooks_{};

    };

    template<typename T>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t>
    CheckpointablePassthrough<T>::registerCheckpointTarget(T0 *target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=](std::string savedStateId) { target->saveState(savedStateId); });
        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() { target->acceptState(); });
        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=](std::string savedStateId) { target->restoreState(savedStateId); });
        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename T>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> CheckpointablePassthrough<T>::registerCacheableCheckpointTarget(T0 *target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=](std::string savedStateId) { target->saveState(savedStateId); });
        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() { target->acceptState(); target->setClean(); });
        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=](std::string savedStateId) { target->restoreState(savedStateId); target->setClean(); });
        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename T>
    void CheckpointablePassthrough<T>::saveState(std::string savedStateId) noexcept {
        if(saved_states_stack_.back().saved_state_id != savedStateId) {
            for (auto &cb : pre_save_hooks_) {
                cb();
            }

            this->underlying().notify_save_state(savedStateId);
            saved_states_stack_.emplace_back(savedStateId);

            for (auto &cb : post_save_hooks_) {
                cb();
            }
        }
    }

    template<typename T>
    void CheckpointablePassthrough<T>::restoreState(std::string savedStateId) noexcept {
        if(saved_states_stack_.back().saved_state_id == savedStateId) {
            for (auto &cb : pre_restore_hooks_) {
                cb();
            }

            this->underlying().notify_restore_state(savedStateId);
            saved_states_stack_.pop_back();

            for (auto &cb : post_restore_hooks_) {
                cb();
            }
        }
    }

    template<typename T>
    void CheckpointablePassthrough<T>::acceptState() noexcept {
        for (auto &cb: pre_accept_hooks_) {
            cb();
        }

        this->underlying().notify_accept_state();
        this->saved_states_stack_.clear();

        for (auto &cb: post_accept_hooks_) {
            cb();
        }
    }

}


#endif //TRANSMISSION_NETWORKS_APP_CHECKPOINTABLEPASSTHROUGH_H
