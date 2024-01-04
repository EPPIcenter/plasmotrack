//
// Created by Maxwell Murphy on 4/24/23.
//


#ifndef TRANSMISSION_NETWORKS_APP_DISCRETEPDF_H
#define TRANSMISSION_NETWORKS_APP_DISCRETEPDF_H

#include "core/computation/PartialLikelihood.h"
#include "core/distributions/DiscreteDistribution.h"
#include "core/parameters/Parameter.h"

#include <cmath>
#include <memory>

namespace transmission_nets::core::distributions {

    template<typename T>
    class DiscretePDF : public computation::PartialLikelihood {
    public:
        using p_Parameter = std::shared_ptr<core::parameters::Parameter<T>>;
        DiscretePDF(p_Parameter target, std::shared_ptr<DiscreteDistribution> probabilities, std::string label = "");
        computation::Likelihood value() override;
        std::string identifier() override;

    private:
        p_Parameter target_;
        std::shared_ptr<DiscreteDistribution> probabilities_;
        std::string label_;
    };

    template<typename T>
    DiscretePDF<T>::DiscretePDF(
            DiscretePDF<T>::p_Parameter target,
            std::shared_ptr<DiscreteDistribution> probabilities,
            std::string label
            ) : target_(std::move(target)), probabilities_(std::move(probabilities)), label_(std::move(label))
    {
        target_->registerCacheableCheckpointTarget(this);
        target_->add_post_change_listener([=, this]() { this->setDirty(); });

        probabilities_->registerCacheableCheckpointTarget(this);
        probabilities_->add_post_change_listener([=, this]() { this->setDirty(); });
    }

    template<typename T>
    computation::Likelihood DiscretePDF<T>::value() {
        if (isDirty()) {
            int idx = std::round(target_->value());

            if (idx < 0 || idx >= (int) probabilities_->value().size()) {
                value_ = -std::numeric_limits<computation::Likelihood>::infinity();
            } else {
                value_ = std::log(probabilities_->value().at(idx));
            }

            setClean();
        }
        return value_;
    }

    template<typename T>
    std::string DiscretePDF<T>::identifier() {
        return {"DiscretePDF"};
    }



}// namespace transmission_nets::core::distributions

#endif//TRANSMISSION_NETWORKS_APP_DISCRETEPDF_H