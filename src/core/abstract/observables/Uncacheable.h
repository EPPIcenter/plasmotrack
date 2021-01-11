//
// Created by Maxwell Murphy on 3/11/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_UNCACHEABLE_H
#define TRANSMISSION_NETWORKS_APP_UNCACHEABLE_H

#include "core/abstract/crtp.h"
#include "core/abstract/observables/Observable.h"


namespace transmission_nets::core::abstract {

    template<typename T, typename ValueType>
    class Uncacheable : public crtp<T, Uncacheable, ValueType> {
        using ChangeCallback = std::function<void()>;
        CRTP_CREATE_EVENT(pre_change, ChangeCallback)
        CRTP_CREATE_EVENT(post_change, ChangeCallback)

    public:

        void setValue(ValueType value) noexcept;

        void initializeValue(ValueType value) noexcept;

        const ValueType& value() const noexcept;

    };

    template<typename T, typename ValueType>
    void Uncacheable<T, ValueType>::setValue(ValueType const value) noexcept {
        assert(this->underlying().isSaved());
        this->underlying().notify_pre_change();
        this->underlying().value_ = value;
        this->underlying().notify_post_change();
    }

    template<typename T, typename ValueType>
    void Uncacheable<T, ValueType>::initializeValue(ValueType const value) noexcept {
        this->underlying().value_ = value;
    }

    template<typename T, typename ValueType>
    [[nodiscard]] const ValueType& Uncacheable<T, ValueType>::value() const noexcept {
        return this->underlying().value_;
    }

}



#endif //TRANSMISSION_NETWORKS_APP_UNCACHEABLE_H
