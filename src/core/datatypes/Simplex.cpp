//
// Created by Maxwell Murphy on 4/17/20.
//

#include "Simplex.h"

Simplex::Simplex(const unsigned int totalElements) : total_elements_(totalElements) {
    coefficients_.resize(total_elements_);
    assert(total_elements_ > 0);
    coefficients_.setOnes();
    coefficients_ = coefficients_ / coefficients_.sum();
}

Simplex::Simplex(const std::initializer_list<double> freqs) : total_elements_(freqs.size()) {
    coefficients_.resize(total_elements_);
    assert(total_elements_ > 0);
    coefficients_.setZero();
    set(freqs);
}

Simplex::Simplex(const std::vector<double> freqs) : total_elements_(freqs.size()) {
    coefficients_.resize(total_elements_);
    assert(total_elements_ > 0);
    coefficients_.setZero();
    set(freqs);
}

Simplex::Simplex(DynamicArray freqs) : total_elements_(freqs.size()), coefficients_(freqs) {}


void Simplex::set(const std::vector<double> valueArray) {
    assert(valueArray.size() == total_elements_);
    for (unsigned int i = 0; i < total_elements_; ++i) {
        coefficients_(i) = valueArray[i];
    }
    coefficients_ = coefficients_ / coefficients_.sum();
}


void Simplex::set(const unsigned int idx, const double value) {
    assert(idx < total_elements_);
//    assert(value < 1);
    coefficients_(idx) = 0;
    coefficients_ = (coefficients_ / coefficients_.sum()) * (1 - value);
    coefficients_(idx) = value;
}

double Simplex::frequencies(const unsigned int idx) const noexcept {
    return coefficients_(idx);
}

const DynamicArray& Simplex::frequencies() const noexcept {
    return coefficients_;
}

unsigned int Simplex::totalElements() const noexcept {
    return total_elements_;
}

std::ostream &operator<<(std::ostream &os, const Simplex &simplex) {
    os << "frequencies: " << simplex.coefficients_;
    return os;
}

double Simplex::min() const noexcept {
    return coefficients_.minCoeff();
}

double Simplex::max() const noexcept {
    return coefficients_.maxCoeff();
}