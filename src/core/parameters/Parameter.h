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

#include "core/utils/forwarding_utils.h"


template<typename T>
class Parameter : public Observable<Parameter<T>>,
                  public Checkpointable<Parameter<T>, T> {

    using ChangeCallback = std::function<void()>;
    CREATE_EVENT(pre_change, ChangeCallback);
    CREATE_EVENT(post_change, ChangeCallback);

public:

    template <typename Args, ENABLE_IF(NonSelf<Args, Parameter<T>>())>
    explicit Parameter(Args&& args) : value_(std::forward<Args>(args)) {
        std::cout << "parameter forwarded c'tor" << std::endl;
    };

    template<typename T0>
    Parameter(const std::initializer_list<T0> il) : value_(il) {
        std::cout << "initializer list constructor" << std::endl;
    };

    Parameter() : value_({}) {
        std::cout << "parameter empty c'tor" << std::endl;
    };

    void setValue(T const value) noexcept {
        assert(this->isSaved());
        this->notify_pre_change();
        value_ = std::move(value);
        this->notify_post_change();
    }

    void setLabel(const std::string& label) noexcept {
        label_ = label;
    }

    [[nodiscard]] std::string label() const noexcept {
        return label_;
    }

    [[nodiscard]] constexpr T value() const noexcept {
        return value_;
    }


protected:
    friend class Checkpointable<Parameter<T>, T>;

    T value_;
    std::string label_{};
};

#endif //TRANSMISSION_NETWORKS_APP_PARAMETER3_H
