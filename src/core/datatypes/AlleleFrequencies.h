//
// Created by Maxwell Murphy on 3/8/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIES_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIES_H

template<int MaxAlleles>
class AlleleFrequenciesVector {
public:

private:
    unsigned int total_alleles_;
    ProbabilityVector<MaxAlleles> allele_frequencies_{};
};

#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCIES_H
