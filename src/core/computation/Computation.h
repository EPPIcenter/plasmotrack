//
// Created by Maxwell Murphy on 1/27/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMPUTATION_H
#define TRANSMISSION_NETWORKS_APP_COMPUTATION_H

namespace transmission_nets::core::computation {

    // A value computed from inputs. value() evaluates it; peek() is the value as last evaluated, which a Checkpointable
    // restores along with the computation, so it is what a running sum over computations holds for this one (see
    // Accumulator). peek() reads nothing and computes nothing.
    template<typename T>
    class Computation {
    public:
        T peek() const noexcept {
            return value_;
        }

        virtual T value() = 0;

        virtual ~Computation() = default;

    protected:
        T value_{};
    };

}// namespace transmission_nets::core::computation

#endif//TRANSMISSION_NETWORKS_APP_COMPUTATION_H
