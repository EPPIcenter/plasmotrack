//
// Created by Maxwell Murphy on 1/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CACHEABLE_H
#define TRANSMISSION_NETWORKS_APP_CACHEABLE_H

#include <functional>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

template<typename T>
class Cacheable : public crtp<T, Cacheable> {
    using SetDirtyCallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(set_dirty, SetDirtyCallbackType)

public:
    bool isDirty() noexcept;

    void setDirty() noexcept;

    void setClean() noexcept;

    template<typename T0>
    void registerDirtyTarget(T0 *target);

protected:
    bool is_dirty_{true};
};

template<typename T>
bool Cacheable<T>::isDirty() noexcept {
    return this->underlying().is_dirty_;
}

template<typename T>
void Cacheable<T>::setDirty() noexcept {
    this->underlying().is_dirty_ = true;
    this->notify_set_dirty();
}

template<typename T>
void Cacheable<T>::setClean() noexcept {
    this->underlying().is_dirty_ = false;
}

template<typename T>
template<typename T0>
void Cacheable<T>::registerDirtyTarget(T0 *target) {
    this->underlying().add_set_dirty_listener([=]() { target->setDirty(); });
}

#endif //TRANSMISSION_NETWORKS_APP_CACHEABLE_H
