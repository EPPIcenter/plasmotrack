//
// Created by Maxwell Murphy on 2/23/20.
//

#include <utility>
#include <string>

#include "Locus.h"

Locus::Locus(std::string label) : uid(newUID++), label(std::move(label)) {}

unsigned int Locus::newUID = 0;

bool Locus::operator<(const Locus &rhs) const {
    return uid < rhs.uid;
}

bool Locus::operator>(const Locus &rhs) const {
    return rhs < *this;
}

bool Locus::operator<=(const Locus &rhs) const {
    return !(rhs < *this);
}

bool Locus::operator>=(const Locus &rhs) const {
    return !(*this < rhs);
}
