//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SIMPLEX_H
#define TRANSMISSION_NETWORKS_APP_SIMPLEX_H

#include <ostream>
#include "core/datatypes/Matrix.h"


// TODO: might want to make this a general frequency data structure, nothing particularly special about alleles

template<int MaxElements>
class Simplex {
public:
    explicit Simplex(unsigned int totalElements);

    Simplex(std::initializer_list<double> freqs);

    explicit Simplex(std::vector<double> freqs);

    void set(std::vector<double> valueArray);

//    void set(unsigned int idx, double value);

//    ProbabilityVector<MaxElements> frequencies() const noexcept {
//        return frequencies_;
//    }

    [[nodiscard]] double frequencies(const unsigned int idx) const noexcept {
        assert(idx < total_elements_);
        return frequencies_(idx);
    };

    [[nodiscard]] unsigned int totalElements() const noexcept {
        return total_elements_;
    };

    friend std::ostream &operator<<(std::ostream &os, const Simplex &vector) {
        for (unsigned int i = 0; i < vector.total_elements_; ++i) {
            os << vector.frequencies_[i] << ", ";
        }
        return os;
    };

private:
    unsigned int total_elements_;
    ProbabilityVector<MaxElements> frequencies_{};
};

template<int MaxAlleles>
Simplex<MaxAlleles>::Simplex(const unsigned int totalElements) : total_elements_(totalElements) {
    frequencies_.setZero();
    frequencies_(0) = 1;
}

template<int MaxAlleles>
Simplex<MaxAlleles>::Simplex(const std::initializer_list<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.setZero();
    set(freqs);
}

template<int MaxAlleles>
Simplex<MaxAlleles>::Simplex(const std::vector<double> freqs) : total_elements_(freqs.size()) {
    frequencies_.setZero();
    set(freqs);
}

template<int MaxAlleles>
void Simplex<MaxAlleles>::set(const std::vector<double> valueArray) {
    assert(valueArray.size() == total_elements_);
    for (unsigned int i = 0; i < total_elements_; ++i) {
        frequencies_[i] = valueArray[i];
    }
    frequencies_ = frequencies_ / frequencies_.sum();
}

//template<int MaxElements>
//void Simplex<MaxElements>::set(unsigned int idx, double value) {
//    assert(idx < total_elements_);
//    assert(value <= 1);
//
//
//}

#endif //TRANSMISSION_NETWORKS_APP_SIMPLEX_H
