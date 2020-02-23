//
// Created by Maxwell Murphy on 1/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CACHEABLE_H
#define TRANSMISSION_NETWORKS_APP_CACHEABLE_H

#include "core/abstract/crtp.h"
#include "Observable.h"

template<typename T>
class Cacheable : public crtp<T, Cacheable> {
    using SetDirtyCallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(set_dirty, SetDirtyCallbackType);

public:
    bool isDirty() noexcept {
        return this->underlying().is_dirty_;
    }

    void setDirty() noexcept {
        this->underlying().is_dirty_ = true;
        this->notify_set_dirty();
    }

    void setClean() noexcept {
        this->underlying().is_dirty_ = false;
    }

protected:
    bool is_dirty_{true};
};

#endif //TRANSMISSION_NETWORKS_APP_CACHEABLE_H
