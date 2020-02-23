//
// Created by Maxwell Murphy on 1/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_DATA_H
#define TRANSMISSION_NETWORKS_APP_DATA_H

#include <utility>

template <typename T>
class Data {

public:
    Data(std::string id, T value) : value_(value), id_(std::move(id)) {};

    template<typename T0>
    Data(std::string id, T0 &&value) : value_(std::forward<T0>(value)), id_(std::move(id)) {};

    [[nodiscard]] constexpr T value() const noexcept {
        return this->value_;
    }

    [[nodiscard]] std::string id() const noexcept {
        return id_;
    }

private:
    T value_;
    std::string id_;
};

#endif //TRANSMISSION_NETWORKS_APP_DATA_H
