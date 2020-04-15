//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H


class AbstractSampler {
public:
    virtual ~AbstractSampler() = default;
    virtual void update() noexcept = 0;
    virtual void adapt() noexcept {};
    virtual void adapt([[maybe_unused]] unsigned int idx) noexcept {};
};


#endif //TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
