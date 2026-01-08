//
// Created by Maxwell Murphy on 1/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H

#include <utility>
#include <vector>
#include <string>
// #include <source_location>
#include <functional>

// #include <boost/container/flat_set.hpp>
// #include <boost/pool/pool_alloc.hpp>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {

//     template <typename T>
//     consteval auto func_name() {
//         const auto& loc = std::source_location::current();
//         return loc.function_name();
//     }
//
//     template <typename T>
//     consteval std::string_view type_of_impl_() {
//         constexpr std::string_view functionName = func_name<T>();
// //        std::string_view sub_type;
// //
// //        for (unsigned int ii = 72; ii < functionName.size(); ++ii) {
// //            if (functionName[ii] == '<') {
// //                sub_type = functionName.substr(72, ii);
// //                break;
// //            }
// //        }
//
// //        return sub_type;
//         return functionName;
//     }
//
//     template <typename T>
//     constexpr auto type_of(T&& arg) {
//         return type_of_impl_<decltype(arg)>();
//     }
//
//     template <typename T>
//     constexpr auto type_of() {
//         return type_of_impl_<T>();
//     }
//
    /*
     * CRTP mixin to enable checkpointing of a value. Allows for the underlying class with field value_ to be saved and restored.
     */
    template<typename T, typename ValueType>
    class Checkpointable : public crtp<T, Checkpointable, ValueType> {
        using CallbackType            = std::function<void()>;
        using AcceptCallbackType      = std::function<void()>;
        using SaveRestoreCallbackType = std::function<void(int saved_state_id)>;
        CRTP_CREATE_EVENT(save_state, SaveRestoreCallbackType)
        CRTP_CREATE_EVENT(accept_state, AcceptCallbackType)
        CRTP_CREATE_EVENT(restore_state, SaveRestoreCallbackType)

        struct StateCheckpoint {
            StateCheckpoint(ValueType savedState, const int savedStateId) : saved_state(savedState), saved_state_id(savedStateId) {};
            ValueType saved_state;
            int saved_state_id;
        };

        // static constexpr std::string_view checkpointable_type = type_of<T>();

    public:
        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0* target);

        template<typename T0>
        std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCacheableCheckpointTarget(T0* target);

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
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCheckpointTarget(T0* target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=, this](int savedStateId) {
            target->saveState(savedStateId);
        });

        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=, this]() {
            target->acceptState();
        });

        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=, this](int savedStateId) {
            target->restoreState(savedStateId);
        });


        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }


    template<typename T, typename ValueType>
    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCacheableCheckpointTarget(T0* target) {
        ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=](int savedStateId) {
            target->saveState(savedStateId);
        });

        ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() {
            target->acceptState();
            target->setClean();
        });

        ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=](int savedStateId) {
            target->restoreState(savedStateId);
            target->setClean();
        });


        return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
    }

    template<typename T, typename ValueType>
    void Checkpointable<T, ValueType>::saveState(int savedStateId) noexcept {
//        fmt::print("Saving state for {} with id {}\n", checkpointable_type, savedStateId);
        if (!isSaved() or saved_states_stack_.back().saved_state_id != savedStateId) {
            for (auto& cb : pre_save_hooks_) {
                cb(savedStateId);
            }

            this->underlying().notify_save_state(savedStateId);

            // fmt::print("Saving state for {} with id {}\n", checkpointable_type, savedStateId);
            // saved_states_stack_.emplace_back(std::move(this->underlying().value()), savedStateId);
            auto val = this->underlying().value();
            auto tmp = StateCheckpoint(val, savedStateId);
            saved_states_stack_.push_back(tmp);
            // fmt::print("Saved state for {} with id {}\n", checkpointable_type, savedStateId);

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
