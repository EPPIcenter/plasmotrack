//
// Created by Maxwell Murphy on 1/27/20.
//

#include "AlleleCounts.h"

AlleleCounts &AlleleCounts::operator-=(const AlleleCounts &rhs) {
    true_positive_count -= rhs.true_positive_count;
    true_negative_count -= rhs.true_negative_count;
    false_positive_count -= rhs.false_positive_count;
    false_negative_count -= rhs.false_negative_count;
    return *this;
}

AlleleCounts &AlleleCounts::operator+=(const AlleleCounts &rhs) {
    true_positive_count += rhs.true_positive_count;
    true_negative_count += rhs.true_negative_count;
    false_positive_count += rhs.false_positive_count;
    false_negative_count += rhs.false_negative_count;
    return *this;
}

AlleleCounts::AlleleCounts() {
    true_positive_count = 0;
    true_negative_count = 0;
    false_positive_count = 0;
    false_negative_count = 0;
}

AlleleCounts operator+(AlleleCounts lhs, const AlleleCounts &rhs) {
    return lhs += rhs;
}

AlleleCounts operator-(AlleleCounts lhs, const AlleleCounts &rhs) {
    return lhs -= rhs;
}

std::ostream &operator<<(std::ostream &out, const AlleleCounts &a) {
    return out  << "TPC: " << a.true_positive_count << "\n"
                << "TNC: " << a.true_negative_count << "\n"
                << "FPC: " << a.false_positive_count << "\n"
                << "FNC: " << a.false_negative_count << "\n";
}
