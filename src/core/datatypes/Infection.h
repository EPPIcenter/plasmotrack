//
// Created by Maxwell Murphy on 1/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_INFECTION_H
#define TRANSMISSION_NETWORKS_APP_INFECTION_H


template <typename GeneticImpl>
struct Infection {
    std::vector<Data<GeneticImpl>> observed_genotype;
    std::vector<Parameter<GeneticImpl>> latent_genotype;
};

#endif //TRANSMISSION_NETWORKS_APP_INFECTION_H
