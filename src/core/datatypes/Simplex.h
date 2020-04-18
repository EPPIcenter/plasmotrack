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

    friend std::ostream &operator<<(std::ostream &os, const Simplex &simplex);

    explicit Simplex(const DynamicArray freqs);

    Simplex(std::initializer_list<double> freqs);

    void set(std::vector<double> valueArray);

    void set(unsigned int idx, double value);

    [[nodiscard]] double frequencies(const unsigned int idx) const noexcept;

    [[nodiscard]] DynamicArray frequencies() const noexcept;

    [[nodiscard]] unsigned int totalElements() const noexcept;

    private:
    unsigned int total_elements_;
    DynamicArray frequencies_{};
};




#endif //TRANSMISSION_NETWORKS_APP_SIMPLEX_H
