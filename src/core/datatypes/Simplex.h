//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEX_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEX_H

#include <ostream>
#include <vector>

#include "core/datatypes/Matrix.h"


class Simplex {
public:
    explicit Simplex(unsigned int totalElements);

    explicit Simplex(std::vector<double> freqs);

    explicit Simplex(const DynamicProbabilityVector freqs);

    Simplex(std::initializer_list<double> freqs);

    void set(std::vector<double> valueArray);

    [[nodiscard]] double frequencies(const unsigned int idx) const noexcept;

    [[nodiscard]] unsigned int totalElements() const noexcept;

    friend std::ostream &operator<<(std::ostream &os, const Simplex &vector) {
        for (unsigned int i = 0; i < vector.total_elements_; ++i) {
            os << vector.frequencies_[i] << ", ";
        }
        return os;
    };

private:
    unsigned int total_elements_;
    DynamicProbabilityVector frequencies_{};
};




#endif //TRANSMISSION_NETWORKS_APP_SIMPLEX_H
