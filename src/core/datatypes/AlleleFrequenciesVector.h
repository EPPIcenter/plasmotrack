//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIESVECTOR_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIESVECTOR_H

#include "core/datatypes/Matrix.h"


// TODO: might want to make this a general frequency data structure, nothing particularly special about alleles

template<int MaxAlleles>
class AlleleFrequenciesVector {
public:
    explicit AlleleFrequenciesVector(const unsigned int totalAlleles);

    AlleleFrequenciesVector(const std::initializer_list<double> freqs);

    explicit AlleleFrequenciesVector(const std::vector<double> freqs);;

    void set(const std::vector<double> valueArray);

    ProbabilityVector<MaxAlleles> alleleFrequencies() const noexcept {
        return allele_frequencies_;
    };

    [[nodiscard]] double alleleFrequencies(int idx) const noexcept {
        return allele_frequencies_(idx);
    };

    [[nodiscard]] unsigned int totalAlleles() const noexcept {
        return total_alleles_;
    }

private:
    unsigned int total_alleles_;
    ProbabilityVector<MaxAlleles> allele_frequencies_{};
};

template<int MaxAlleles>
AlleleFrequenciesVector<MaxAlleles>::AlleleFrequenciesVector(const unsigned int totalAlleles) : total_alleles_(totalAlleles) {
    allele_frequencies_.setZero();
    allele_frequencies_(0) = 1;
}

template<int MaxAlleles>
AlleleFrequenciesVector<MaxAlleles>::AlleleFrequenciesVector(const std::initializer_list<double> freqs) : total_alleles_(freqs.size()) {
    allele_frequencies_.setZero();
    set(freqs);
}

template<int MaxAlleles>
AlleleFrequenciesVector<MaxAlleles>::AlleleFrequenciesVector(const std::vector<double> freqs) : total_alleles_(freqs.size()) {
    allele_frequencies_.setZero();
    set(freqs);
}

template<int MaxAlleles>
void AlleleFrequenciesVector<MaxAlleles>::set(const std::vector<double> valueArray) {
    assert(valueArray.size() == total_alleles_);
    for (unsigned int i = 0; i < total_alleles_; ++i) {
        allele_frequencies_[i] = valueArray[i];
    }
    allele_frequencies_ = allele_frequencies_ / allele_frequencies_.sum();
}

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIESVECTOR_H
