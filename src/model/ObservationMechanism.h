//
// Created by Maxwell Murphy on 1/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H

#include "core/datatypes/Data.h"

template <
        typename AlleleStorage,
        template<typename...> typename SavedStatePolicy = std::optional,
        template<typename...> typename DependencyStore = boost::container::flat_set>
class ObservationMechanism : public PartialLogLikelihood<SavedStatePolicy, DependencyStore> {
private:
    Parameter<AlleleStorage>* latent_;
    Data<AlleleStorage>* observed_;
    Parameter<float>* epsilon_pos;
    Parameter<float>* epsilon_neg;

public:
    ObservationMechanism(Parameter<AlleleStorage> *latent, Data<AlleleStorage> *observed) : latent_(latent), observed_(observed) {
        latent->addDependency(this);
    }

    using PartialLogLikelihood<SavedStatePolicy, DependencyStore>::setDirty;

    void setDirty(AbstractComputation *ptr) noexcept override {
        std::cout << "SetDirty" << "\n";
        this->setDirty();
    }

    void setDirty(AbstractParameter *ptr) noexcept override {
        std::cout << "SetDirty" << "\n";
        this->setDirty();
    }

    float value() override {
        return AlleleStorage::falsePositiveCount(latent_->value(), observed_->value());
    };
};

#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISM_H
