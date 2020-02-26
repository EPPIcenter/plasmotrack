//
// Created by Maxwell Murphy on 1/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"


template<typename T, typename ValueType>
class Checkpointable : public crtp<T, Checkpointable, ValueType> {
    using CallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(save_state, CallbackType);
    CRTP_CREATE_EVENT(accept_state, CallbackType);
    CRTP_CREATE_EVENT(restore_state, CallbackType);

public:

    template<typename T0>
    void registerCheckpointTarget(T0& target) {
        this->underlying().add_save_state_listener([&]() { target.saveState();});
        this->underlying().add_accept_state_listener([&]() { target.acceptState();});
        this->underlying().add_restore_state_listener([&]() { target.restoreState();});
    }

    void saveState() noexcept {
        if(!this->saved_state_) {
            std::cout << "Saving State" << std::endl;
            this->underlying().notify_save_state();
            this->saved_state_ = this->underlying().value();
        }
    }

    void restoreState() noexcept {
        if (this->saved_state_) {
            std::cout << "Restoring State" << std::endl;
            this->underlying().notify_restore_state();
            this->underlying().value_ = *(this->saved_state_);
            this->saved_state_.reset();
        }
    }

    void acceptState() noexcept {
        if (this->saved_state_) {
            std::cout << "Accepting State" << std::endl;
            this->underlying().notify_accept_state();
            this->saved_state_.reset();
        }
    }

    bool constexpr isSaved() noexcept {
        return this->saved_state_.has_value();
    }

protected:
    std::optional<ValueType> saved_state_{};
};

#endif //TRANSMISSION_NETWORKS_APP_CHECKPOINTABLE_H
