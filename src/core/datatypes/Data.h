//
// Created by Maxwell Murphy on 1/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_DATA_H
#define TRANSMISSION_NETWORKS_APP_DATA_H

#include <utility>

#include "core/utils/forwarding_utils.h"

namespace transmission_nets::core::datatypes {

    template <typename T>
    class Data {

    public:
        template <typename Args, ENABLE_IF(core::utils::NonSelf<Args, Data<T>>())>
        explicit Data(Args&& args) : value_(std::forward<Args>(args)) {
//        std::cout << "data forward c'tor" << std::endl;
        }

        void setLabel(const std::string& label) noexcept {
            label_ = label;
        }

        [[nodiscard]] T value() const noexcept {
            return value_;
        }

        [[nodiscard]] std::string label() const noexcept {
            return label_;
        }

    private:
        T value_;
        std::string label_;
    };
}





#endif //TRANSMISSION_NETWORKS_APP_DATA_H
