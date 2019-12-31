//
// Created by Maxwell Murphy on 10/16/19.
//

#ifndef TRANSMISSION_NETWORK_MODEL_NODE_H
#define TRANSMISSION_NETWORK_MODEL_NODE_H

#include <bitset>
#include <vector>
#include "TransmissionNetworkModelConfig.h"

typedef std::bitset<MAX_ALLELES> AlleleSet;
typedef std::vector<AlleleSet> Genotype;
typedef std::vector<float> AlleleFrequencies;
typedef std::vector<AlleleFrequencies> LocusAlleleFrequencies;

struct Node {
public:
    Genotype latent_genotype;
};


struct UnobservedNode : Node { };


struct ObservedNode : Node {
public:
    const std::string label;
    const std::vector<bool> has_data;
    const Genotype observed_genotype;
//    const UnobservedNode unobserved_parent;
};

struct SourceNode : Node {
public:
    const std::string label;
    LocusAlleleFrequencies locus_allele_frequencies;
};



#endif //TRANSMISSION_NETWORK_MODEL_NODE_H
