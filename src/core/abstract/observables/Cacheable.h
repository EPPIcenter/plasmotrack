//
// Cacheable: a CRTP mixin that gives a computation a dirty flag and announces each time it is marked dirty, so that
// what depends on it is marked dirty too (see Dependency.h). value() recomputes when dirty and marks clean.
//

#ifndef TRANSMISSION_NETWORKS_APP_CACHEABLE_H
#define TRANSMISSION_NETWORKS_APP_CACHEABLE_H

#include <functional>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {
    template<typename T>
    class Cacheable : public crtp<T, Cacheable> {
        using SetDirtyCallbackType = std::function<void()>;
        CRTP_CREATE_EVENT(set_dirty, SetDirtyCallbackType)

    public:
        bool isDirty() noexcept;

        void setDirty() noexcept;

        void setClean() noexcept;

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
}// namespace transmission_nets::core::abstract


#endif//TRANSMISSION_NETWORKS_APP_CACHEABLE_H
