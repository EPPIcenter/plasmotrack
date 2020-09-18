//
// Created by Maxwell Murphy on 5/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PRIOR_H
#define TRANSMISSION_NETWORKS_APP_PRIOR_H


#include "core/computation/PartialLikelihood.h"


namespace transmission_nets::core::priors {
    template<typename Distribution, typename TargetParam, typename ...Args>
    class Prior : public computation::PartialLikelihood {

    public:
        explicit Prior(TargetParam &target, Args&&... args);

        computation::Likelihood value() override {
            if (isDirty()) {
                this->value_ = log(pdf(dist_, target_.value()));
                this->value_ = std::isnan(this->value_) ? -std::numeric_limits<double>::infinity() : this->value_;
                this->setClean();
            }
            return this->value_;
        }

    private:
        TargetParam &target_;
        Distribution dist_;
    };

    template<typename Distribution, typename TargetParam, typename...Args>
    Prior<Distribution, TargetParam, Args...>::Prior(TargetParam &target, Args&&... args) : target_(target), dist_(std::forward<Args>(args)...) {
        target_.registerCacheableCheckpointTarget(this);
        target_.add_post_change_listener([=, this]() { this->setDirty(); });
        this->setDirty();
    }
}



#endif //TRANSMISSION_NETWORKS_APP_PRIOR_H
