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
        explicit Simplex(unsigned int totalElements);
        Simplex();

        explicit Simplex(const std::vector<double>& freqs);

        Simplex(const std::initializer_list<double>& freqs);

        friend std::ostream& operator<<(std::ostream& os, const Simplex& simplex);

        void set(const std::vector<double>& valueArray);

        void set(unsigned int idx, double value);

        [[nodiscard]] double frequencies(unsigned int idx) const noexcept;

        [[nodiscard]] const std::vector<double>& frequencies() const noexcept;

        [[nodiscard]] unsigned int totalElements() const noexcept;

        [[nodiscard]] double min() const noexcept;

        [[nodiscard]] double max() const noexcept;

        [[nodiscard]] std::string serialize() const noexcept;

    private:
        std::vector<double> coefficients_{};
        unsigned int total_elements_;
        double min_;
        double max_;
    };

}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_SIMPLEX_H
