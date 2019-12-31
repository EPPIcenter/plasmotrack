//
// Created by Maxwell Murphy on 11/6/19.
//

#include "NoSuperInfection.h"

//LogLik NoSuperInfection::calc_llik_transmission(
//        const ObservedNode& parent_node, const ObservedNode& child_node, int num_transmissions,
//        const std::vector<int>& total_alleles, int total_loci) {
//
//    assert(total_loci > 0);
//    assert(num_transmissions > 0);
//
//    // Calculate on a per locus basis the likelihood of the number of mutations and retained alleles given some number
//    // of transmissions
//
//    LogLik llik = 0.0;
//    unsigned int mutations;
//    unsigned int non_mutations;
//    unsigned int parent_coi;
//    unsigned int child_coi;
//    for (int i = 0; i < total_loci; ++i) {
//        if (parent_node.has_data[i] && child_node.has_data[i]) {
//            parent_coi = parent_node.latent_genotype[i].count(); // 1X
//            child_coi = (parent_node.latent_genotype[i] & child_node.latent_genotype[i]).count(); // 11
//            mutations = (~parent_node.latent_genotype[i] & child_node.latent_genotype[i]).count(); // 01
//            non_mutations = (~parent_node.latent_genotype[i] & ~child_node.latent_genotype[i]).count() -
//                            (MAX_ALLELES - total_alleles[i]); // 00
//            llik += log(calc_prob_coi_transmission(parent_coi, child_coi, num_transmissions));
//            llik += mutations * log(transmission_mech_parameters.mutation_rate);
//            llik += non_mutations * log(1 - transmission_mech_parameters.mutation_rate);
//        }
//    }
//    return llik;
//}

COITransitionMatrix NoSuperInfection::zt_multiplicative_binomial(const NoSuperInfectionParameters &params) {
    // Calculate the probability of all possible transitions according to a zero-truncated multiplicative binomial
    // distribution.
    auto tmp = combinations_matrix.array() *
               Eigen::pow(params.transmission_complexity_prob, matrix_one.array()) *
               Eigen::pow(1 - params.transmission_complexity_prob, matrix_two.array()) *
               Eigen::pow(params.transmission_complexity_assoc, matrix_three.array());

    // Normalize the rows so that it is a valid transition matrix
    auto res = tmp.array().colwise() / tmp.rowwise().sum();
    return res;
}

void NoSuperInfection::initialize_matrices() {
    // Initialize matrices for use in calculating the zero-truncated multiplicative binomial probabilities.
    combinations_matrix = TransitionMatrix::Zero();
    matrix_one = TransitionMatrix::Zero();
    matrix_two = TransitionMatrix::Zero();
    matrix_three = TransitionMatrix::Zero();
    for (long f = 0; f < combinations_matrix.rows(); ++f) {
        for (long t = 0; t <= f; ++t) {
            combinations_matrix(f, t) = boost::math::binomial_coefficient<double>(f + 1, t + 1);
            matrix_one(f, t) = t + 1;
            matrix_two(f, t) = f - t;
            matrix_three(f, t) = (t + 1) * (f - t);
        }
    }
}

void NoSuperInfection::set_transmission_mech_parameters(const NoSuperInfectionParameters &parameters) {
    if (parameters.mutation_rate != transmission_mech_parameters.mutation_rate) {
        update_mutation_rate(parameters);
    }

    if (
            parameters.transmission_complexity_prob != transmission_mech_parameters.transmission_complexity_prob ||
            parameters.transmission_complexity_assoc != transmission_mech_parameters.transmission_complexity_assoc
            ) {
        update_coi_transmission_probabilities(parameters);
    }
}

Probability NoSuperInfection::calc_prob_coi_transmission(int from_coi, int to_coi, int num_transmissions) {
    assert(from_coi <= MAX_COI);
    assert(to_coi <= MAX_COI);
    assert(num_transmissions <= MAX_TRANSMISSIONS);

    return transmission_probability[num_transmissions - 1](from_coi - 1, to_coi - 1);
}

void NoSuperInfection::update_mutation_rate(const NoSuperInfectionParameters &parameters) {
    transmission_mech_parameters.mutation_rate = parameters.mutation_rate;
}

void NoSuperInfection::update_coi_transmission_probabilities(const NoSuperInfectionParameters &parameters) {
    transmission_mech_parameters.transmission_complexity_assoc = parameters.transmission_complexity_assoc;
    transmission_mech_parameters.transmission_complexity_prob = parameters.transmission_complexity_prob;

    transmission_probability[0] = zt_multiplicative_binomial(transmission_mech_parameters);
    for (int g = 1; g < MAX_TRANSMISSIONS; ++g) {
        transmission_probability[g] = transmission_probability[0] * transmission_probability[g - 1];
    }
}

LogLik NoSuperInfection::calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set,
                                                int total_loci, const std::vector<int> &total_alleles) {

//    assert(total_loci > 0);
//    assert(num_transmissions > 0);

    // Calculate on a per locus basis the likelihood of the number of mutations and retained alleles given some number
    // of transmissions

    LogLik llik = 0.0;
//    unsigned int mutations;
//    unsigned int non_mutations;
//    unsigned int parent_coi;
//    unsigned int child_coi;
//    for (int i = 0; i < total_loci; ++i) {
//        if (parent_node.has_data[i] && child_node.has_data[i]) {
//            parent_coi = parent_node.latent_genotype[i].count(); // 1X
//            child_coi = (parent_node.latent_genotype[i] & child_node.latent_genotype[i]).count(); // 11
//            mutations = (~parent_node.latent_genotype[i] & child_node.latent_genotype[i]).count(); // 01
//            non_mutations = (~parent_node.latent_genotype[i] & ~child_node.latent_genotype[i]).count() -
//                            (MAX_ALLELES - total_alleles[i]); // 00
//            llik += log(calc_prob_coi_transmission(parent_coi, child_coi, num_transmissions));
//            llik += mutations * log(transmission_mech_parameters.mutation_rate);
//            llik += non_mutations * log(1 - transmission_mech_parameters.mutation_rate);
//        }
//    }
    return llik;
}

LogLik NoSuperInfection::calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set,
                                                const SourceNode &source_node, int total_loci,
                                                const std::vector<int> &total_alleles) {
    return 0;
}

LogLik NoSuperInfection::calc_llik_transmission(const ObservedNode &node, const SourceNode &source_node, int total_loci,
                                                const std::vector<int> &total_alleles) {
    return 0;
}

NoSuperInfection::NoSuperInfection(const NoSuperInfectionParameters &transmission_mech_parameters)
        : transmission_mech_parameters(transmission_mech_parameters) {
    initialize_matrices();
    update_coi_transmission_probabilities(transmission_mech_parameters);
}

NoSuperInfection::NoSuperInfection() {
    initialize_matrices();
}

NoSuperInfectionParameters NoSuperInfection::get_transmission_mech_parameters() {
    return transmission_mech_parameters;
}

