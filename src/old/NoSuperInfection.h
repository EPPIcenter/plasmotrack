//
// Created by Maxwell Murphy on 11/6/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H
#define TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H

#include "TransmissionMechanism.h"

struct NoSuperInfectionParameters {
    Probability mutation_rate;
    Probability transmission_complexity_prob;
    float transmission_complexity_assoc;
};

class NoSuperInfection : public TransmissionMechanism<NoSuperInfectionParameters> {
private:
    NoSuperInfection();

    NoSuperInfection(const NoSuperInfectionParameters &transmission_mech_parameters);

    NoSuperInfectionParameters transmission_mech_parameters;

    TransmissionProbability transmission_probability;
    TransitionMatrix combinations_matrix;
    TransitionMatrix matrix_one;
    TransitionMatrix matrix_two;
    TransitionMatrix matrix_three;

    void initialize_matrices();

    COITransitionMatrix zt_multiplicative_binomial(const NoSuperInfectionParameters &params);

    void update_coi_transmission_probabilities(const NoSuperInfectionParameters &params);

    void update_mutation_rate(const NoSuperInfectionParameters &parameters);

    Probability calc_prob_coi_transmission(int from_coi, int to_coi, int num_transmissions);

public:
    LogLik
    calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set, int total_loci,
                           const std::vector<int> &total_alleles) override;

    LogLik
    calc_llik_transmission(const ObservedNode &node, const std::vector<Node> &parent_set, const SourceNode &source_node,
                           int total_loci, const std::vector<int> &total_alleles) override;

    LogLik calc_llik_transmission(const ObservedNode &node, const SourceNode &source_node, int total_loci,
                                  const std::vector<int> &total_alleles) override;

    void set_transmission_mech_parameters(const NoSuperInfectionParameters &params) override;

    NoSuperInfectionParameters get_transmission_mech_parameters() override;

};



#endif //TRANSMISSION_NETWORKS_APP_NOSUPERINFECTION_H
