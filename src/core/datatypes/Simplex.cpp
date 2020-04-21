//
// Created by Maxwell Murphy on 4/17/20.
//

#include "Simplex.h"

Simplex::Simplex(const unsigned int totalElements) : total_elements_(totalElements) {
    frequencies_.resize(total_elements_);
    assert(total_elements_ > 0);
    frequencies_.setOnes();
    frequencies_ = frequencies_ / frequencies_.sum();
}

Simplex::Simplex(const std::initializer_list<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.resize(total_elements_);
    assert(total_elements_ > 0);
    frequencies_.setZero();
    set(freqs);
}

Simplex::Simplex(const std::vector<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.resize(total_elements_);
    assert(total_elements_ > 0);
    frequencies_.setZero();
    set(freqs);
}

Simplex::Simplex(DynamicArray freqs) : total_elements_(freqs.size()), frequencies_(freqs) {}


void Simplex::set(const std::vector<double> valueArray) {
    assert(valueArray.size() == total_elements_);
    for (unsigned int i = 0; i < total_elements_; ++i) {
        frequencies_(i) = valueArray[i];
    }
    frequencies_ = frequencies_ / frequencies_.sum();
}


void Simplex::set(const unsigned int idx, const double value) {
    assert(idx < total_elements_);
    assert(value < 1);
    frequencies_(idx) = 0;
    frequencies_ = (frequencies_ / frequencies_.sum()) * (1 - value);
    frequencies_(idx) = value;
}

double Simplex::frequencies(const unsigned int idx) const noexcept {
    return frequencies_(idx);
}

DynamicArray Simplex::frequencies() const noexcept {
    return frequencies_;
}

unsigned int Simplex::totalElements() const noexcept {
    return total_elements_;
}

std::ostream &operator<<(std::ostream &os, const Simplex &simplex) {
    os << "frequencies: " << simplex.frequencies_;
    return os;
}

