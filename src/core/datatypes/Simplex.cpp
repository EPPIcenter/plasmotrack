//
// Created by Maxwell Murphy on 4/17/20.
//

#include "Simplex.h"

Simplex::Simplex(const unsigned int totalElements) : total_elements_(totalElements) {
    frequencies_.resize(total_elements_);
    frequencies_.setZero();
    frequencies_(0) = 1;
}

Simplex::Simplex(const std::initializer_list<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.resize(total_elements_);
    frequencies_.setZero();
    set(freqs);
}

Simplex::Simplex(const std::vector<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.resize(total_elements_);
    frequencies_.setZero();
    set(freqs);
}

void Simplex::set(const std::vector<double> valueArray) {
    assert(valueArray.size() == total_elements_);
    for (unsigned int i = 0; i < total_elements_; ++i) {
        frequencies_(i) = valueArray[i];
    }
    frequencies_ = frequencies_ / frequencies_.sum();
}

double Simplex::frequencies(const unsigned int idx) const noexcept {
    return frequencies_(idx);
}

unsigned int Simplex::totalElements() const noexcept {
    return total_elements_;
}

Simplex::Simplex(DynamicProbabilityVector freqs) : total_elements_(freqs.size()), frequencies_(freqs) {}
