//
// Created by Maxwell Murphy on 12/6/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV2_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV2_H

namespace transmission_nets::model::observation_process {

    /*
     * Observation process to incorporate error.
     * This implementation will scale the observation error values by the total number of alleles, so
     * instead of a probability, it's a rate of observation error. This makes the value insensitive to the
     * total number of alleles present at a marker.
     */

    template<typename GeneticsImpl>
    class ObservationProcessLikelihoodv2 : public core::computation::PartialLikelihood {
        static constexpr auto truePositiveCount  = &GeneticsImpl::truePositiveCount;
        static constexpr auto falsePositiveCount = &GeneticsImpl::falsePositiveCount;
        static constexpr auto trueNegativeCount  = &GeneticsImpl::trueNegativeCount;
        static constexpr auto falseNegativeCount = &GeneticsImpl::falseNegativeCount;

    public:
        using p_Parameterdouble = std::shared_ptr<core::parameters::Parameter<double>>;
        ObservationProcessLikelihoodv2(
                std::shared_ptr<core::datatypes::Data<GeneticsImpl>> observedGenetics,
                std::shared_ptr<core::parameters::Parameter<GeneticsImpl>> latentGenetics,
                p_Parameterdouble expectedFalsePositives,
                p_Parameterdouble expectedFalseNegatives);

        core::computation::Likelihood value() override;
        std::string identifier() override;

    private:
        std::shared_ptr<core::datatypes::Data<GeneticsImpl>> observed_genetics_;
        std::shared_ptr<core::parameters::Parameter<GeneticsImpl>> latent_genetics_;
        p_Parameterdouble expected_false_positives_;
        p_Parameterdouble expected_false_negatives_;
        unsigned int total_alleles_;
    };

    template<typename GeneticsImpl>
    ObservationProcessLikelihoodv2<GeneticsImpl>::ObservationProcessLikelihoodv2(
            std::shared_ptr<core::datatypes::Data<GeneticsImpl>> observedGenetics,
            std::shared_ptr<core::parameters::Parameter<GeneticsImpl>> latentGenetics,
            p_Parameterdouble expectedFalsePositives,
            p_Parameterdouble expectedFalseNegatives) : observed_genetics_(std::move(observedGenetics)),
                                                        latent_genetics_(std::move(latentGenetics)),
                                                        expected_false_positives_(std::move(expectedFalsePositives)),
                                                        expected_false_negatives_(std::move(expectedFalseNegatives)) {
        total_alleles_ = latent_genetics_->value().totalAlleles();
        latent_genetics_->add_post_change_listener([=, this]() { this->setDirty(); });
        latent_genetics_->registerCacheableCheckpointTarget(this);

        expected_false_positives_->add_post_change_listener([=, this]() { this->setDirty(); });
        expected_false_positives_->registerCacheableCheckpointTarget(this);

        expected_false_negatives_->add_post_change_listener([=, this]() { this->setDirty(); });
        expected_false_negatives_->registerCacheableCheckpointTarget(this);

        this->setDirty();
        this->value();
    }

    template<typename GeneticsImpl>
    std::string ObservationProcessLikelihoodv2<GeneticsImpl>::identifier() {
        return {"ObservationProcessLikelihoodv2"};
    }

    template<typename GeneticsImpl>
    core::computation::Likelihood ObservationProcessLikelihoodv2<GeneticsImpl>::value() {
        if (this->isDirty()) {
            const auto& latent_genetics = latent_genetics_->value();
            const auto& observed_genetics = observed_genetics_->value();
            const int true_positive_count  = truePositiveCount(latent_genetics, observed_genetics);
            const int false_positive_count = falsePositiveCount(latent_genetics, observed_genetics);
            const int true_negative_count  = trueNegativeCount(latent_genetics, observed_genetics);
            const int false_negative_count = falseNegativeCount(latent_genetics, observed_genetics);

            const double expected_false_positives = expected_false_positives_->value();
            const double expected_false_negatives = expected_false_negatives_->value();

            value_ = true_positive_count * log(1 - (expected_false_positives / total_alleles_)) +
                     true_negative_count * log(1 - (expected_false_negatives / total_alleles_)) +
                     false_positive_count * log(expected_false_positives / total_alleles_) +
                     false_negative_count * log(expected_false_negatives / total_alleles_);

            this->setClean();
        }

        assert(value_ < std::numeric_limits<core::computation::Likelihood>::infinity());

        return value_;
    }


}// namespace transmission_nets::model::observation_process


#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONPROCESSLIKELIHOODV2_H
