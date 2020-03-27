//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H


class AbstractSampler {
public:
    virtual ~AbstractSampler() = default;
    virtual void update() = 0;

};


#endif //TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
