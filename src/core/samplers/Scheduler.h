//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SCHEDULER_H

#include <vector>

#include "AbstractSampler.h"

class Scheduler {
public:
    void registerSampler(AbstractSampler* sampler) {
        samplers.push_back(sampler);
    }

    void updateAndAdapt() {
        for(const auto sampler : samplers) {
            sampler->update();
        }
        for(const auto sampler : samplers) {
            sampler->adapt();
        }
    }

    void update() {
        for(const auto sampler : samplers) {
            sampler->update();
        }
    }

private:

    std::vector<AbstractSampler*> samplers{};
    std::vector<double> weights{};
};


#endif //TRANSMISSION_NETWORKS_APP_SCHEDULER_H
