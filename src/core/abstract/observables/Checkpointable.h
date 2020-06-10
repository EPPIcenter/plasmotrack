//
// Created by Maxwell Murphy on 1/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H

#include <type_traits>
#include <optional>


#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

template<typename T, typename ValueType>
class Checkpointable : public crtp<T, Checkpointable, ValueType> {
    using CallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(save_state, CallbackType)
    CRTP_CREATE_EVENT(accept_state, CallbackType)
    CRTP_CREATE_EVENT(restore_state, CallbackType)

public:

    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCheckpointTarget(T0* target);

    template<typename T0>
    std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> registerCacheableCheckpointTarget(T0* target);

    void saveState() noexcept;

    void restoreState() noexcept;

    void acceptState() noexcept;

    bool constexpr isSaved() noexcept;

protected:
    std::optional<ValueType> saved_state_{};
};

template<typename T, typename ValueType>
template<typename T0>
std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCheckpointTarget(T0 *target) {
    ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=]() { target->saveState(); });
    ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() { target->acceptState(); });
    ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=]() { target->restoreState(); });
    return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
}

template<typename T, typename ValueType>
template<typename T0>
std::tuple<ListenerId_t, ListenerId_t, ListenerId_t> Checkpointable<T, ValueType>::registerCacheableCheckpointTarget(T0 *target) {
    ListenerId_t saveStateEventId = this->underlying().add_save_state_listener([=]() {
        target->saveState();
    });

    ListenerId_t acceptStateEventId = this->underlying().add_accept_state_listener([=]() {
        target->acceptState();
        target->value();
        target->setClean();
    });

    ListenerId_t restoreStateEventId = this->underlying().add_restore_state_listener([=]() {
        target->restoreState();
        target->value();
        target->setClean();
    });

    return std::make_tuple(saveStateEventId, acceptStateEventId, restoreStateEventId);
}

template<typename T, typename ValueType>
void Checkpointable<T, ValueType>::saveState() noexcept {
    if (!this->saved_state_) {
        this->underlying().notify_save_state();
        this->saved_state_ = this->underlying().value();
    }
}

template<typename T, typename ValueType>
void Checkpointable<T, ValueType>::restoreState() noexcept {
    if (this->saved_state_) {
        this->underlying().notify_restore_state();
        this->underlying().value_ = *(this->saved_state_);
        this->saved_state_.reset();
    }
}

template<typename T, typename ValueType>
void Checkpointable<T, ValueType>::acceptState() noexcept {
    if (this->saved_state_) {
        this->underlying().notify_accept_state();
        this->saved_state_.reset();
    }
}

template<typename T, typename ValueType>
constexpr bool Checkpointable<T, ValueType>::isSaved() noexcept {
    return this->saved_state_.has_value();
}

#endif //TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
