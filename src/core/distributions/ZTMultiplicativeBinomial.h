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

#include "core/datatypes/Matrix.h"

#include "core/parameters/Parameter.h"

#include "core/computation/Computation.h"


namespace transmission_nets::core::distributions {
    template<int MAX_COUNT>
    class ZTMultiplicativeBinomial : public computation::Computation<datatypes::TransitionMatrix<MAX_COUNT + 1>>,
                                     public abstract::Observable<ZTMultiplicativeBinomial<MAX_COUNT>>,
                                     public abstract::Cacheable<ZTMultiplicativeBinomial<MAX_COUNT>>,
                                     public abstract::Checkpointable<ZTMultiplicativeBinomial<MAX_COUNT>, datatypes::TransitionMatrix<MAX_COUNT + 1>> {
    public:
        ZTMultiplicativeBinomial(parameters::Parameter<double> &prob,
                                 parameters::Parameter<double> &assoc);

        datatypes::TransitionMatrix<MAX_COUNT + 1> value() noexcept override;


    private:
        friend class abstract::Checkpointable<ZTMultiplicativeBinomial<MAX_COUNT>, datatypes::TransitionMatrix<MAX_COUNT + 1>>;
        friend class abstract::Cacheable<ZTMultiplicativeBinomial<MAX_COUNT>>;

        static const datatypes::SquareMatrix<double, MAX_COUNT + 1> mat1;
        static const datatypes::SquareMatrix<double, MAX_COUNT + 1> mat2;
        static const datatypes::SquareMatrix<double, MAX_COUNT + 1> mat3;
        static const datatypes::SquareMatrix<double, MAX_COUNT + 1> combo_matrix;

        parameters::Parameter<double> &prob_;
        parameters::Parameter<double> &assoc_;

    };


    template<int MAX_COUNT>
    datatypes::SquareMatrix<double, MAX_COUNT + 1> initMat1() {
        datatypes::SquareMatrix<double, MAX_COUNT + 1> a;
        for (int j = 0; j <= MAX_COUNT; ++j) {
            for (int k = 0; k <= j; ++k) {
                a(j, k) = k;
            }
        }
        return a;
    }

    template<int MAX_COUNT>
    datatypes::SquareMatrix<double, MAX_COUNT + 1> initMat2() {
        datatypes::SquareMatrix<double, MAX_COUNT + 1> a;
        for (int j = 0; j <= MAX_COUNT; ++j) {
            for (int k = 0; k <= j; ++k) {
                a(j, k) = j - k;
            }
        }
        return a;
    }

    template<int MAX_COUNT>
    datatypes::SquareMatrix<double, MAX_COUNT + 1> initMat3() {
        datatypes::SquareMatrix<double, MAX_COUNT + 1> a;
        for (int j = 0; j <= MAX_COUNT; ++j) {
            for (int k = 0; k <= j; ++k) {
                a(j, k) = (k) * (j - k);
            }
        }
        return a;
    }

    template<int MAX_COUNT>
    datatypes::SquareMatrix<double, MAX_COUNT + 1> initCombosMat() {
        datatypes::SquareMatrix<double, MAX_COUNT + 1> a;
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
    const datatypes::SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat1 = initMat1<MAX_COUNT>();


    template<int MAX_COUNT>
    const datatypes::SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat2 = initMat2<MAX_COUNT>();


    template<int MAX_COUNT>
    const datatypes::SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::mat3 = initMat3<MAX_COUNT>();


    template<int MAX_COUNT>
    const datatypes::SquareMatrix<double, MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::combo_matrix = initCombosMat<MAX_COUNT>();


    template<int MAX_COUNT>
    ZTMultiplicativeBinomial<MAX_COUNT>::ZTMultiplicativeBinomial(parameters::Parameter<double> &prob, parameters::Parameter<double> &assoc): prob_(prob), assoc_(assoc) {
        prob_.registerCacheableCheckpointTarget(this);
        prob_.add_post_change_listener([=, this]() { this->setDirty(); });

        assoc_.registerCacheableCheckpointTarget(this);
        assoc_.add_post_change_listener([=, this]() { this->setDirty(); });
        this->setDirty();
        this->value();
    }


    template<int MAX_COUNT>
    datatypes::TransitionMatrix<MAX_COUNT + 1> ZTMultiplicativeBinomial<MAX_COUNT>::value() noexcept {
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
}


#endif //TRANSMISSION_NETWORKS_APP_ZTMULTIPLICATIVEBINOMIAL_H
