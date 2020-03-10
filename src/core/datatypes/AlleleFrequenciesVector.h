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
    explicit AlleleFrequenciesVector(const unsigned int totalAlleles) : total_alleles_(totalAlleles) {
        allele_frequencies_.setZero();
        allele_frequencies_(0) = 1;
    }

    AlleleFrequenciesVector(const std::initializer_list<double> freqs) : total_alleles_(freqs.size()) {
        allele_frequencies_.setZero();
        set(freqs);
    };

    explicit AlleleFrequenciesVector(const std::vector<double> freqs) : total_alleles_(freqs.size()) {
        allele_frequencies_.setZero();
        set(freqs);
    };

    //
//    AlleleFrequenciesVector(const AlleleFrequenciesVector<MaxAlleles> &other) : total_alleles_(
//            std::move(other.total_alleles_)), allele_frequencies_(std::move(other.allele_frequencies_)) {
//        std::cout << "allele frequencies move constructor" << std::endl;
//    }

    void set(const std::vector<double> valueArray) {
        assert(valueArray.size() == total_alleles_);
        for (unsigned int i = 0; i < total_alleles_; ++i) {
            allele_frequencies_[i] = valueArray[i];
        }
        allele_frequencies_ = allele_frequencies_ / allele_frequencies_.sum();
    }

    ProbabilityVector<MaxAlleles> alleleFrequencies() const noexcept {
        return allele_frequencies_;
    };

    [[nodiscard]] double alleleFrequencies(int idx) const noexcept {
        return allele_frequencies_(idx);
    };

private:
    unsigned int total_alleles_;
    ProbabilityVector<MaxAlleles> allele_frequencies_{};
};

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIESVECTOR_H
