//
// Created by Maxwell Murphy on 12/9/19.
//

#include <iostream>

#include "Likelihood.h"


float Likelihood::value() {
    if (!is_dirty_) {
        return value_;
    }

    for (const auto& lik: dirty_dependencies_) {
        value_ -= lik->peek();
        value_ += lik->value();
    }

    setClean();

    return value_;
}


void Likelihood::pushLikelihood(PartialLikelihood* ptr){
    likelihoods_.insert(ptr);
    dirty_dependencies_.insert(ptr);
    value_ += ptr->value();
    ptr->addDependency(this);
}

void Likelihood::reinitialize() {
    value_ = 0;
    for (const auto& lik: likelihoods_) {
        value_ += lik->value();
    }
    setClean();
}

Likelihood::Likelihood() {
    value_ = 0;
}

float Likelihood::peek() {
    return value_;
}
