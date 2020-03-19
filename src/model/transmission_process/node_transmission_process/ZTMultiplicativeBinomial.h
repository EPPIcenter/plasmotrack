//
// Created by Maxwell Murphy on 2/7/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ZTMULTIPLICATIVEBINOMIAL_H
#define TRANSMISSION_NETWORKS_APP_ZTMULTIPLICATIVEBINOMIAL_H

#include <boost/math/special_functions/binomial.hpp>
#include <Eigen/Core>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/datatypes/TransitionMatrix.h"

#include "core/parameters/Parameter.h"

#include "core/computation/Computation.h"


template<int MAX_COUNT>
class ZTMultiplicativeBinomial : public Computation<ProbabilityMatrix <MAX_COUNT + 1>>,
                                 public Observable<ZTMultiplicativeBinomial<MAX_COUNT>>,
                                 public Cacheable<ZTMultiplicativeBinomial<MAX_COUNT>>,
                                 public Checkpointable<ZTMultiplicativeBinomial<MAX_COUNT>, ProbabilityMatrix <MAX_COUNT + 1>> {
public:
    ZTMultiplicativeBinomial(Parameter<double> &prob,
                             Parameter<double> &assoc);

    ProbabilityMatrix <MAX_COUNT + 1> value() noexcept override;


private:
    friend class Checkpointable<ZTMultiplicativeBinomial<MAX_COUNT>, ProbabilityMatrix <MAX_COUNT + 1>>;
    friend class Cacheable<ZTMultiplicativeBinomial<MAX_COUNT>>;

    static const SquareMatrix<double, MAX_COUNT + 1> mat1;
    static const SquareMatrix<double, MAX_COUNT + 1> mat2;
    static const SquareMatrix<double, MAX_COUNT + 1> mat3;
    static const SquareMatrix<double, MAX_COUNT + 1> combo_matrix;

    Parameter<double> &prob_;
    Parameter<double> &assoc_;

};


template<int MAX_COUNT>
SquareMatrix<double, MAX_COUNT + 1> initMat1() {
    SquareMatrix<double, MAX_COUNT + 1> a;
    for (int j = 0; j <= MAX_COUNT; ++j) {
        for (int k = 0; k <= j; ++k) {
            a(j, k) = k;
        }
    }
    return a;
}

template<int MAX_COUNT>
SquareMatrix<double, MAX_COUNT + 1> initMat2() {
    SquareMatrix<double, MAX_COUNT + 1> a;
    for (int j = 0; j <= MAX_COUNT; ++j) {
        for (int k = 0; k <= j; ++k) {
            a(j, k) = j - k;
        }
    }
    return a;
}

template<int MAX_COUNT>
SquareMatrix<double, MAX_COUNT + 1> initMat3() {
    SquareMatrix<double, MAX_COUNT + 1> a;
    for (int j = 0; j <= MAX_COUNT; ++j) {
        for (int k = 0; k <= j; ++k) {
            a(j, k) = (k) * (j - k);
        }
    }
    return a;
}

template<int MAX_COUNT>
SquareMatrix<double, MAX_COUNT + 1> initCombosMat() {
    SquareMatrix<double, MAX_COUNT + 1> a;
    for (int j = 0; j <= MAX_COUNT; ++j) {
        for (int k = 0; k <= j; ++k) {
            if(!j or !k) {
                a(j, k) = 0;
                continue;
            }
            a(j, k) = boost::math::binomial_coefficient<double>(j, k);
        }
    }
    return a;
}


template<int MAX_COUNT>
const SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat1 = initMat1<MAX_COUNT>();


template<int MAX_COUNT>
const SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat2 = initMat2<MAX_COUNT>();


template<int MAX_COUNT>
const SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat3 = initMat3<MAX_COUNT>();


template<int MAX_COUNT>
const SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::combo_matrix = initCombosMat<MAX_COUNT>();


template<int MAX_COUNT>
ZTMultiplicativeBinomial<MAX_COUNT>::ZTMultiplicativeBinomial(Parameter<double> &prob, Parameter<double> &assoc): prob_(prob), assoc_(assoc) {
    prob_.registerCacheableCheckpointTarget(*this);
    prob_.add_post_change_listener([&]() { this->setDirty(); });
    
    assoc_.registerCacheableCheckpointTarget(*this);
    assoc_.add_post_change_listener([&]() { this->setDirty(); });
    }


template<int MAX_COUNT>
ProbabilityMatrix <MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::value() noexcept {
    if (this->isDirty()) {
        auto tmp = combo_matrix.array() *
                   Eigen::pow(prob_.value(), mat1.array()) *
                   Eigen::pow(1 - prob_.value(), mat2.array()) *
                   Eigen::pow(assoc_.value(), mat3.array());
        this->value_ = tmp.array().colwise() / tmp.rowwise().sum();
        this->value_.row(0).setZero();
        this->setClean();
    }
    return this->value_;
}

#endif //TRANSMISSION_NETWORKS_APP_ZTMULTIPLICATIVEBINOMIAL_H
