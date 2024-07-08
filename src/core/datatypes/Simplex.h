//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEX_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEX_H


#include "core/datatypes/Matrix.h"


#include <fmt/core.h>

#include <vector>


namespace transmission_nets::core::datatypes {

    class Simplex {
    public:
        explicit Simplex(unsigned char totalElements);
        Simplex();

        explicit Simplex(const std::vector<double>& freqs);

        Simplex(const std::initializer_list<double>& freqs);

        friend std::ostream& operator<<(std::ostream& os, const Simplex& simplex);

        void set(const std::vector<double>& valueArray);

        void set(unsigned char idx, double value);

        [[nodiscard]] double frequencies(unsigned char idx) const noexcept;

        [[nodiscard]] std::vector<double> frequencies() const noexcept;

        [[nodiscard]] unsigned char totalElements() const noexcept;

        [[nodiscard]] double min() const noexcept;

        [[nodiscard]] double max() const noexcept;

        [[nodiscard]] std::string serialize() const noexcept;

    private:
        // std::vector<float> coefficients_{};
        std::array<double, 64> coefficients_{};
        unsigned char total_elements_;
        double min_;
        double max_;
    };

}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_SIMPLEX_H
