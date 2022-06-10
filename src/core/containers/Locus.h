//
// Created by Maxwell Murphy on 2/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOCUS_H
#define TRANSMISSION_NETWORKS_APP_LOCUS_H

#include <string>

namespace transmission_nets::core::containers {

    class Locus {
    public:
        const unsigned int uid;
        const std::string label;

        explicit Locus(std::string label, int total_alleles);

        virtual ~Locus();

        [[nodiscard]] unsigned int totalAlleles() const noexcept;

        bool operator<(const Locus& rhs) const noexcept;

        bool operator>(const Locus& rhs) const noexcept;

        bool operator<=(const Locus& rhs) const noexcept;

        bool operator>=(const Locus& rhs) const noexcept;

    private:
        static unsigned int newUID;
        unsigned int total_alleles_;
    };

}// namespace transmission_nets::core::containers


#endif//TRANSMISSION_NETWORKS_APP_LOCUS_H
