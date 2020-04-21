//
// Created by Maxwell Murphy on 1/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARAMETER3_H
#define TRANSMISSION_NETWORKS_APP_PARAMETER3_H

#include <boost/container/flat_map.hpp>
#include <functional>
#include <optional>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Uncacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/utils/forwarding_utils.h"


template<typename T>
class Parameter : public Observable<Parameter<T>>,
                  public Uncacheable<Parameter<T>, T>,
                  public Checkpointable<Parameter<T>, T> {


public:

    template <typename Args, ENABLE_IF(NonSelf<Args, Parameter<T>>())>
    explicit Parameter(Args&& args) : value_(std::forward<Args>(args)) {
        std::cout << "parameter forwarded c'tor" << std::endl;
    };

    template<typename T0>
    Parameter(const std::initializer_list<T0> il) : value_(il) {
        std::cout << "initializer list constructor" << std::endl;
    };

    Parameter() {
        std::cout << "parameter empty c'tor" << std::endl;
    };


    void setLabel(const std::string& label) noexcept {
        label_ = label;
    }

    [[nodiscard]] std::string label() const noexcept {
        return label_;
    }


protected:
    friend class Checkpointable<Parameter<T>, T>;
    friend class Uncacheable<Parameter<T>, T>;

    T value_;
    std::string label_{};
};

#endif //TRANSMISSION_NETWORKS_APP_PARAMETER3_H
