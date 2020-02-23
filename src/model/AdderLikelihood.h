//
// Created by Maxwell Murphy on 1/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H

#include "../core/PartialLogLikelihood.h"
#include "../core/Parameter3.h"

template <template<typename...> typename SavedStatePolicy = std::optional,
        template<typename...> typename DependencyStore = boost::container::flat_set>
class AdderLikelihood : public PartialLogLikelihood<SavedStatePolicy, DependencyStore> {
private:
    Parameter<float>* a_;
    Parameter<float>* b_;

public:
    AdderLikelihood(Parameter<float> *a, Parameter<float> *b) : PartialLogLikelihood<SavedStatePolicy, DependencyStore>(), a_(a), b_(b) {
        a_->addDependency(this);
        b_->addDependency(this);
    }

    using PartialLogLikelihood<SavedStatePolicy, DependencyStore>::setDirty;

    void setDirty([[maybe_unused]] AbstractComputation *ptr) noexcept override {
        this->setDirty();
    }

    void setDirty(AbstractParameter *ptr) noexcept override {
        if(ptr == a_) {
            std::cout << "A updated" << "\n";
        } else if (ptr == b_) {
            std::cout << "B updated" << "\n";
        }

        this->setDirty();
    }

    [[nodiscard]] float value() noexcept override {
        if(this->isDirty()) {
            this->value_ = a_->value() + b_->value();
            this->setClean();
        }
        return this->value_;
    }
};


#endif //TRANSMISSION_NETWORKS_APP_ADDERLIKELIHOOD_H
