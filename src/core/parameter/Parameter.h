//
// Created by Maxwell Murphy on 1/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARAMETER3_H
#define TRANSMISSION_NETWORKS_APP_PARAMETER3_H

#include <boost/container/flat_map.hpp>
#include <functional>
#include <optional>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Checkpointable.h"


template<typename T>
class Parameter : public Observable<Parameter<T>>,
                  public Checkpointable<Parameter<T>, T> {

    using ChangeCallback = std::function<void()>;
    CREATE_EVENT(pre_change, ChangeCallback);
    CREATE_EVENT(post_change, ChangeCallback);

public:

    Parameter(std::string id, T value) : value_(std::move(value)), id_(std::move(id)) {};

    void setValue(T const value) noexcept {
        assert(this->isSaved());
        this->notify_pre_change();
        value_ = std::move(value);
        this->notify_post_change();
    }

    [[nodiscard]] constexpr T value() const noexcept {
        return value_;
    }


protected:
    friend class Checkpointable<Parameter<T>, T>;

    T value_;
    const std::string id_;
};

#endif //TRANSMISSION_NETWORKS_APP_PARAMETER3_H
