//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELECOUNTS_H
#define TRANSMISSION_NETWORKS_APP_ALLELECOUNTS_H

#include <iostream>

namespace transmission_nets::model::observation_process {

    struct AlleleCounts {
        unsigned int true_positive_count{};
        unsigned int true_negative_count{};
        unsigned int false_positive_count{};
        unsigned int false_negative_count{};

        AlleleCounts();

        AlleleCounts& operator+=(const AlleleCounts& rhs);

        AlleleCounts& operator-=(const AlleleCounts& rhs);
    };

    AlleleCounts operator+(AlleleCounts lhs, const AlleleCounts& rhs);

    AlleleCounts operator-(AlleleCounts lhs, const AlleleCounts& rhs);

    std::ostream& operator<<(std::ostream& out, const AlleleCounts& a);

}// namespace transmission_nets::model::observation_process


#endif//TRANSMISSION_NETWORKS_APP_ALLELECOUNTS_H
