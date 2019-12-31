//
// Created by Maxwell Murphy on 12/10/19.
//

#include "AdderLikelihood.h"

float AdderLikelihood::value() {
    if(!isDirty()) {
        return value_;
    }

    value_ = a_->value() + b_->value();
    setClean();

    return value_;
}

AdderLikelihood::AdderLikelihood(IntegerParameter *a, IntegerParameter *b)
        {
    a_ = a;
    b_ = b;
    a_->addDependency(this);
    b_->addDependency(this);
}

float AdderLikelihood::peek() {
    return value_;
}
