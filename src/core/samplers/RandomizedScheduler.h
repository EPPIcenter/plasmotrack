//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H

#include <boost/random.hpp>

#include "Scheduler.h"



template<typename Engine=boost::random::mt19937>
class RandomizedScheduler : public Scheduler {

public:
    explicit RandomizedScheduler(Engine *rng);

    void update() override;

private:
    Engine* rng_;

};


template<typename Engine>
RandomizedScheduler<Engine>::RandomizedScheduler(Engine *rng):rng_(rng) {}

template<typename Engine>
void RandomizedScheduler<Engine>::update() {
    auto indices = randomSequence(0, samplers_.size(), rng_);
    for (const auto idx : indices) {
        samplers_.at(idx)->update();
    }
}


#endif //TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
