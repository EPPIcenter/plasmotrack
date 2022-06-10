//
// Created by Maxwell Murphy on 2/23/20.
//

#include <string>
#include <utility>

#include "Locus.h"


namespace transmission_nets::core::containers {
    Locus::Locus(std::string label, int total_alleles) : uid(newUID++), label(std::move(label)),
                                                         total_alleles_(total_alleles) {}
    Locus::~Locus() = default;

    unsigned int Locus::newUID = 0;

    unsigned int Locus::totalAlleles() const noexcept {
        return total_alleles_;
    }

    bool Locus::operator<(const Locus& rhs) const noexcept {
        return uid < rhs.uid;
    }

    bool Locus::operator>(const Locus& rhs) const noexcept {
        return rhs < *this;
    }

    bool Locus::operator<=(const Locus& rhs) const noexcept {
        return !(rhs < *this);
    }

    bool Locus::operator>=(const Locus& rhs) const noexcept {
        return !(*this < rhs);
    }


}// namespace transmission_nets::core::containers
