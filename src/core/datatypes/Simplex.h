//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEX_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEX_H

#include <vector>

#include "core/datatypes/Matrix.h"

namespace transmission_nets::core::datatypes {

    class Simplex {
    public:
        explicit Simplex(unsigned int totalElements);

        explicit Simplex(std::vector<double> freqs);

        explicit Simplex(DynamicArray freqs);

        Simplex(std::initializer_list<double> freqs);

        friend std::ostream &operator<<(std::ostream &os, const Simplex &simplex);

        void set(std::vector<double> valueArray);

        void set(unsigned int idx, double value);

        [[nodiscard]] double frequencies(unsigned int idx) const noexcept;

        [[nodiscard]] const DynamicArray& frequencies() const noexcept;

        [[nodiscard]] unsigned int totalElements() const noexcept;

        [[nodiscard]] double min() const noexcept;

        [[nodiscard]] double max() const noexcept;

        [[nodiscard]] const std::string serialize() const noexcept;

    private:
        unsigned int total_elements_;
        DynamicArray coefficients_{};
        Eigen::IOFormat fmt{4, 0, ",", ","};
    };

}




#endif //TRANSMISSION_NETWORKS_APP_SIMPLEX_H
