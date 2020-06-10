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
        samplers_.push_back(sampler);
    }

    virtual void updateAndAdapt() {
        update();
        adapt();
    }

    virtual void update() {
        for(const auto& sampler : samplers_) {
            sampler->update();
        }
    }

    virtual void adapt() {
        for(const auto& sampler : samplers_) {
            sampler->adapt();
        }
    }

    [[nodiscard]] const std::vector<AbstractSampler*>& samplers() const noexcept {
        return samplers_;
    }

protected:

    std::vector<AbstractSampler*> samplers_{};
    std::vector<double> weights_{};
};


#endif //TRANSMISSION_NETWORKS_APP_SCHEDULER_H
