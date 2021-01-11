//
// Created by Maxwell Murphy on 3/10/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CACHEABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_CACHEABLEPASSTHROUGH_H

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

namespace transmission_nets::core::abstract {

    template<typename T>
    class CacheablePassthrough : public crtp<T, CacheablePassthrough> {
        using SetDirtyCallbackType = std::function<void()>;
    CRTP_CREATE_EVENT(set_dirty, SetDirtyCallbackType);

    public:
        void setDirty() noexcept;

        template<typename T0>
        void registerDirtyTarget(T0& target);

    };

    template<typename T>
    void CacheablePassthrough<T>::setDirty() noexcept {
        this->notify_set_dirty();
    }

    template<typename T>
    template<typename T0>
    void CacheablePassthrough<T>::registerDirtyTarget(T0 &target) {
        this->underlying().add_set_dirty_listener([=, this]() { target.setDirty();});
    }

}




#endif //TRANSMISSION_NETWORKS_APP_CACHEABLEPASSTHROUGH_H
