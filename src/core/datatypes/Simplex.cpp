//
// Created by Maxwell Murphy on 4/17/20.
//

#include "Simplex.h"

#include "core/io/serialize.h"


namespace transmission_nets::core::datatypes {

    Simplex::Simplex(const unsigned int totalElements) : total_elements_(totalElements) {
        assert(total_elements_ > 0);
        coefficients_.assign(total_elements_, 1.0 / total_elements_);
        min_ = coefficients_[0];
        max_ = coefficients_[0];
    }

    Simplex::Simplex(const std::initializer_list<float>& freqs) : total_elements_(freqs.size()), min_(std::numeric_limits<float>::max()), max_(std::numeric_limits<float>::min()) {
        coefficients_.resize(total_elements_);
        assert(total_elements_ > 0);
        set(freqs);
    }

    Simplex::Simplex(const std::vector<float>& freqs) : total_elements_(freqs.size()), min_(std::numeric_limits<float>::max()), max_(std::numeric_limits<float>::min()) {
        coefficients_.resize(total_elements_);
        assert(total_elements_ > 0);
        set(freqs);
    }

    void Simplex::set(const std::vector<float>& valueArray) {
        assert(valueArray.size() == total_elements_);
        min_       = std::numeric_limits<float>::max();
        max_       = std::numeric_limits<float>::min();
        float sum = 0;
        for (unsigned int ii = 0; ii < total_elements_; ++ii) {
            coefficients_[ii] = valueArray[ii];
            sum += valueArray[ii];
            min_ = std::min(min_, coefficients_[ii]);
            max_ = std::max(max_, coefficients_[ii]);
        }

        if (sum != 1.0) {
            min_ = std::numeric_limits<float>::max();
            max_ = std::numeric_limits<float>::min();
            for (unsigned int ii = 0; ii < total_elements_; ++ii) {
                coefficients_[ii] = coefficients_[ii] / sum;
                min_              = std::min(min_, coefficients_[ii]);
                max_              = std::max(max_, coefficients_[ii]);
            }
        }
    }

    void Simplex::set(const unsigned int idx, const float value) {
        assert(idx < total_elements_);
        min_               = std::numeric_limits<float>::max();
        max_               = std::numeric_limits<float>::min();
        float prev_value  = coefficients_[idx];
        coefficients_[idx] = 0;
        for (unsigned int ii = 0; ii < total_elements_; ++ii) {
            coefficients_[ii] = (coefficients_[ii] / (1 - prev_value)) * (1 - value);
            min_              = std::min(min_, coefficients_[ii]);
            max_              = std::max(max_, coefficients_[ii]);
        }
        coefficients_[idx] = value;
        min_               = std::min(min_, coefficients_[idx]);
        max_               = std::max(max_, coefficients_[idx]);
    }

    float Simplex::frequencies(const unsigned int idx) const noexcept {
        return coefficients_[idx];
    }

    const std::vector<float>& Simplex::frequencies() const noexcept {
        return coefficients_;
    }

    unsigned int Simplex::totalElements() const noexcept {
        return total_elements_;
    }

    std::ostream& operator<<(std::ostream& os, const Simplex& simplex) {
        os << "frequencies: " << core::io::serialize(simplex.coefficients_);
        return os;
    }

    float Simplex::min() const noexcept {
        return min_;
    }

    float Simplex::max() const noexcept {
        return max_;
    }

    std::string Simplex::serialize() const noexcept {
        std::stringstream ss;
        ss << core::io::serialize(coefficients_);
        return ss.str();
    }

}// namespace transmission_nets::core::datatypes
