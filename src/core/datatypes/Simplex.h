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

        explicit Simplex(const std::vector<float>& freqs);

        Simplex(const std::initializer_list<float>& freqs);

        friend std::ostream& operator<<(std::ostream& os, const Simplex& simplex);

        void set(const std::vector<float>& valueArray);

        void set(unsigned int idx, float value);

        [[nodiscard]] float frequencies(unsigned int idx) const noexcept;

        [[nodiscard]] const std::vector<float>& frequencies() const noexcept;

        [[nodiscard]] unsigned int totalElements() const noexcept;

        [[nodiscard]] float min() const noexcept;

        [[nodiscard]] float max() const noexcept;

        [[nodiscard]] std::string serialize() const noexcept;

    private:
        std::vector<float> coefficients_{};
        unsigned int total_elements_;
        float min_;
        float max_;
    };

}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_SIMPLEX_H
