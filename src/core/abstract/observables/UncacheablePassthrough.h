//
// Created by Maxwell Murphy on 3/11/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_UNCACHEABLEPASSTHROUGH_H
#define TRANSMISSION_NETWORKS_APP_UNCACHEABLEPASSTHROUGH_H

#include <functional>

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"

template<typename T>
class UncacheablePassthrough : public crtp<T, UncacheablePassthrough> {
    using ChangeCallback = std::function<void()>;
    CRTP_CREATE_EVENT(pre_change, ChangeCallback)
    CRTP_CREATE_EVENT(post_change, ChangeCallback)
//    virtual
};

#endif //TRANSMISSION_NETWORKS_APP_UNCACHEABLEPASSTHROUGH_H
