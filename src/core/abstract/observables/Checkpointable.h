//
// Created by Maxwell Murphy on 1/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H

#include <utility>
#include <vector>

#include <boost/pool/pool_alloc.hpp>
#include <boost/container/flat_set.hpp>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {

    /*
     * CRTP mixin to enable checkpointing of a value. Allows for the underlying class with field value_ to be saved and restored.
     */

    template<typename T, typename ValueType>
    class Checkpointable : public crtp<T, Checkpointable, ValueType> {
        using CallbackType = std::function<void()>;
        using AcceptCallbackType = std::function<void()>;
        using SaveRestoreCallbackType = std::function<void(const std::string& saved_state_id)>;
        CRTP_CREATE_EVENT(save_state, SaveRestoreCallbackType)
        CRTP_CREATE_EVENT(accept_state, AcceptCallbackType)
        CRTP_CREATE_EVENT(restore_state, SaveRestoreCallbackType)

        struct StateCheckpoint {
            StateCheckpoint(ValueType savedState, std::string savedStateId) : saved_state(savedState), saved_state_id(std::move(savedStateId)) {}
            ValueType saved_state;
            std::string saved_state_id;
        };

    public:
        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0* target);

        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCacheableCheckpointTarget(T0* target);

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

        void saveState(const std::string& savedStateId) noexcept;

        void restoreState(const std::string& savedStateId) noexcept;

        void acceptState() noexcept;

        bool constexpr isSaved() noexcept;

    protected:
        std::vector<StateCheckpoint, boost::pool_allocator<StateCheckpoint>> saved_states_stack_{};
        std::vector<CallbackType> pre_save_hooks_{};
        std::vector<CallbackType> post_save_hooks_{};
        std::vector<CallbackType> pre_restore_hooks_{};
        std::vector<CallbackType> post_restore_hooks_{};
        std::vector<CallbackType> pre_accept_hooks_{};
        std::vector<CallbackType> post_accept_hooks_{};
    };


    template<typename T, typename ValueType>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCheckpointTarget(T0 *target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=, this](const std::string& savedStateId) {
          target->saveState(savedStateId);
        });

        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=, this]() {
          target->acceptState();
        });

        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=, this](const std::string& savedStateId) {
          target->restoreState(savedStateId);
        });


        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }


    template<typename T, typename ValueType>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCacheableCheckpointTarget(T0 *target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=](const std::string& savedStateId) {
          target->saveState(savedStateId);
        });

        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() {
          target->acceptState();
          target->setClean();
        });

        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=](const std::string& savedStateId) {
          target->restoreState(savedStateId);
          target->setClean();
        });


        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::saveState(const std::string& savedStateId) noexcept {
        if(!isSaved() or saved_states_stack_.back().saved_state_id != savedStateId) {
            for (auto &cb: pre_save_hooks_) {
                cb();
            }

            this->underlying().notify_save_state(savedStateId);
            saved_states_stack_.emplace_back(this->underlying().value(), savedStateId);

            for (auto &cb: post_save_hooks_) {
                cb();
            }
        }
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::restoreState(const std::string& savedStateId) noexcept {
        if(isSaved() and saved_states_stack_.back().saved_state_id == savedStateId) {
            for (auto &cb: pre_restore_hooks_) {
                cb();
            }

            this->underlying().notify_restore_state(savedStateId);
            this->underlying().value_ = std::move(saved_states_stack_.back().saved_state);
            saved_states_stack_.pop_back();

            for (auto &cb: post_restore_hooks_) {
                cb();
            }
        }
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::acceptState() noexcept {
        if (isSaved()) {
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

    template<typename T, typename ValueType>
    constexpr bool Checkpointable<T, ValueType>::isSaved() noexcept {
        return !(this->saved_states_stack_.empty());
    }

}



#endif //TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
