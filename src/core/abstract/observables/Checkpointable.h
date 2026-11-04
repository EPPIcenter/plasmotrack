//
// Checkpointable: a CRTP mixin that gives a value a stack of saved states, so that a Proposal's change to it can be
// accepted or restored, and announces each save, accept and restore to whatever depends on it (see Dependency.h).
// A save is ignored while the same id is on top of the stack, a restore pops only if its id is on top, and an accept
// clears the stack; Proposal.h gives ids and the order of calls.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace transmission_nets::core::abstract {

    template<typename T, typename ValueType>
    class Checkpointable : public crtp<T, Checkpointable, ValueType> {
        using CallbackType            = std::function<void()>;
        using AcceptCallbackType      = std::function<void()>;
        using SaveRestoreCallbackType = std::function<void(int saved_state_id)>;
        CRTP_CREATE_EVENT(save_state, SaveRestoreCallbackType)
        CRTP_CREATE_EVENT(accept_state, AcceptCallbackType)
        CRTP_CREATE_EVENT(restore_state, SaveRestoreCallbackType)

        struct StateCheckpoint {
            StateCheckpoint(ValueType savedState, const int savedStateId) : saved_state(std::move(savedState)), saved_state_id(savedStateId) {};
            ValueType saved_state;
            int saved_state_id;
        };

    public:
        void addPreSaveHook(const SaveRestoreCallbackType& cb) {
            this->pre_save_hooks_.emplace_back(cb);
        }

        void addPostSaveHook(const SaveRestoreCallbackType& cb) {
            this->post_save_hooks_.emplace_back(cb);
        }

        void addPreRestoreHook(const SaveRestoreCallbackType& cb) {
            this->pre_restore_hooks_.emplace_back(cb);
        }

        void addPostRestoreHook(const SaveRestoreCallbackType& cb) {
            this->post_restore_hooks_.emplace_back(cb);
        }

        void addPreAcceptHook(const AcceptCallbackType& cb) {
            this->pre_accept_hooks_.emplace_back(cb);
        }

        void addPostAcceptHook(const AcceptCallbackType& cb) {
            this->post_accept_hooks_.emplace_back(cb);
        }

        // Pushes the current value under `savedStateId`, evaluating it first, unless that id is already on top.
        void saveState(int savedStateId) noexcept;

        void restoreState(int savedStateId) noexcept;

        void acceptState() noexcept;

        bool constexpr isSaved() noexcept;

    protected:
        std::vector<StateCheckpoint> saved_states_stack_{};
        std::vector<SaveRestoreCallbackType> pre_save_hooks_{};
        std::vector<SaveRestoreCallbackType> post_save_hooks_{};
        std::vector<SaveRestoreCallbackType> pre_restore_hooks_{};
        std::vector<SaveRestoreCallbackType> post_restore_hooks_{};
        std::vector<AcceptCallbackType> pre_accept_hooks_{};
        std::vector<AcceptCallbackType> post_accept_hooks_{};
    };


    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::saveState(int savedStateId) noexcept {
        if (!isSaved() or saved_states_stack_.back().saved_state_id != savedStateId) {
            for (auto& cb : pre_save_hooks_) {
                cb(savedStateId);
            }

            this->underlying().notify_save_state(savedStateId);

            // value() returns a copy, which is moved onto the stack.
            saved_states_stack_.emplace_back(this->underlying().value(), savedStateId);

            for (auto& cb : post_save_hooks_) {
                cb(savedStateId);
            }
        }
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::restoreState(int savedStateId) noexcept {
        if (isSaved() and saved_states_stack_.back().saved_state_id == savedStateId) {
            for (auto& cb : pre_restore_hooks_) {
                cb(savedStateId);
            }

            this->underlying().notify_restore_state(savedStateId);
            this->underlying().value_ = std::move(saved_states_stack_.back().saved_state);
            saved_states_stack_.pop_back();

            for (auto& cb : post_restore_hooks_) {
                cb(savedStateId);
            }
        }
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::acceptState() noexcept {
        if (isSaved()) {
            for (auto& cb : pre_accept_hooks_) {
                cb();
            }

            this->underlying().notify_accept_state();
            this->saved_states_stack_.clear();

            for (auto& cb : post_accept_hooks_) {
                cb();
            }
        }
    }

    template<typename T, typename ValueType>
    constexpr bool Checkpointable<T, ValueType>::isSaved() noexcept {
        return !(this->saved_states_stack_.empty());
    }

}// namespace transmission_nets::core::abstract


#endif//TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
