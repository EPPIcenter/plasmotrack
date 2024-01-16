//
// Created by mmurphy on 10/22/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_TEMPERED_H
#define TRANSMISSION_NETWORKS_APP_TEMPERED_H


#include <fmt/core.h>

#include "core/computation/Computation.h"
#include "core/computation/PartialLikelihood.h"


namespace transmission_nets::core::computation {
    template<typename Input>
    struct Tempered : public PartialLikelihood {

        std::shared_ptr<Input> target_;
        float temperature_;

        explicit Tempered(std::shared_ptr<Input> target, float temperature) : target_(std::move(target)), temperature_(temperature) {
            target_->add_set_dirty_listener([=, this]() {
                this->setDirty();
            });
            target_->registerCacheableCheckpointTarget(this);
        }

        Likelihood value() override {
            if (this->isDirty()) {
                this->value_ = temperature_ * target_->value();
                this->setClean();
            }
            return this->value_;
        }

        void setTemperature(float temperature) {
            temperature_ = temperature;
            this->setDirty();
        }

        [[nodiscard]] float getTemperature() const {
            return temperature_;
        }

        std::string identifier() override {
            return fmt::format("Tempered({}, {})", target_->identifier(), temperature_);
        }

    };
}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_TEMPERED_H
