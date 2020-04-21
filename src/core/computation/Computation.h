//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMPUTATION_H
#define TRANSMISSION_NETWORKS_APP_COMPUTATION_H

template <typename T>
class Computation {
public:
    [[nodiscard]] T peek() noexcept {
        return value_;
    };

    virtual T value() = 0;

    virtual ~Computation() = default;

protected:
    T value_{};
};

#endif //TRANSMISSION_NETWORKS_APP_COMPUTATION_H
