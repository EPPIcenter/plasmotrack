//
// Created by Maxwell Murphy on 1/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARAMETER_H
#define TRANSMISSION_NETWORKS_APP_PARAMETER_H

#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Uncacheable.h"

#include "core/utils/forwarding_utils.h"


namespace transmission_nets::core::parameters {

    template<typename T>
    class Parameter : public abstract::Observable<Parameter<T>>,
                      public abstract::Uncacheable<Parameter<T>, T>,
                      public abstract::Checkpointable<Parameter<T>, T> {

    public:
        template<typename Args, ENABLE_IF(core::utils::NonSelf<Args, Parameter<T>>())>
        explicit Parameter(Args&& args) : value_(std::forward<Args>(args)) {}

        template<typename T0>
        Parameter(const std::initializer_list<T0> il) : value_(il) {}

        Parameter() = default;


        void setLabel(const std::string& label) noexcept {
            label_ = label;
        }

        [[nodiscard]] std::string label() const noexcept {
            return label_;
        }

    protected:
        friend class abstract::Checkpointable<Parameter<T>, T>;
        friend class abstract::Uncacheable<Parameter<T>, T>;

        T value_{};
        std::string label_{};
    };

}// namespace transmission_nets::core::parameters


#endif//TRANSMISSION_NETWORKS_APP_PARAMETER_H
