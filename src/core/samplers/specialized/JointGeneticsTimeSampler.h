//
// Created by Maxwell Murphy on 11/8/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_JOINTGENETICSTIMESAMPLER_H
#define TRANSMISSION_NETWORKS_APP_JOINTGENETICSTIMESAMPLER_H


#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/genetics/ParentSetGenotypeProposal.h"
#include "core/samplers/genetics/ParentSetProposalProbability.h"
#include "core/utils/numerics.h"

#include <memory>

namespace transmission_nets::core::samplers::specialized {

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    class JointGeneticsTimeSampler : public AbstractSampler {
    public:
        JointGeneticsTimeSampler(
                std::shared_ptr<InfectionEventImpl> infection,
                std::shared_ptr<ParentSetImpl> ps,
                std::shared_ptr<InfectionEventImpl> latent_parent,
                std::shared_ptr<parameters::Parameter<double>> infection_time,
                std::shared_ptr<T> target,
                std::shared_ptr<Engine> rng,
                double lower_bound = 1.0,
                double upper_bound = 1000.0,
                double variance = 3.0,
                double fn_rate = 0.01,
                double fp_rate = 0.01
                ) noexcept;

        void update() noexcept override;

        [[nodiscard]] unsigned int acceptances() noexcept;

        [[nodiscard]] unsigned int rejections() noexcept;

        [[nodiscard]] double acceptanceRate() noexcept;

        [[nodiscard]] Likelihood calculateSamplingProb() noexcept;

    private:
        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<ParentSetImpl> ps_;
        std::shared_ptr<InfectionEventImpl> latent_parent_;
        std::shared_ptr<parameters::Parameter<double>> infection_duration_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double lower_bound_;
        double upper_bound_;
        double variance_;
        double fn_rate_;
        double fp_rate_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::normal_distribution<> normal_dist_{0, 1};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;

        genetics::ParentSetProposalProbability<InfectionEventImpl, MaxParentSetSize> proposalProbability_{fn_rate_, fp_rate_};
    };

    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    JointGeneticsTimeSampler<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::JointGeneticsTimeSampler(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<ParentSetImpl> ps,
            std::shared_ptr<InfectionEventImpl> latent_parent,
            std::shared_ptr<parameters::Parameter<double>> infection_time,
            std::shared_ptr<T> target,
            std::shared_ptr<Engine> rng,
            double lower_bound,
            double upper_bound,
            double variance,
            double fn_rate,
            double fp_rate) noexcept :
                                                    infection_(std::move(infection)),
                                                    ps_(std::move(ps)),
                                                    latent_parent_(std::move(latent_parent)),
                                                    infection_duration_(std::move(infection_time)),
                                                    target_(std::move(target)),
                                                    rng_(rng),
                                                    lower_bound_(lower_bound),
                                                    upper_bound_(upper_bound),
                                                    variance_(variance),
                                                    fn_rate_(fn_rate),
                                                    fp_rate_(fp_rate)
    {}


    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    void JointGeneticsTimeSampler<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::update() noexcept {
        SAMPLER_STATE_ID stateId = SAMPLER_STATE_ID::JointGeneticsTimeID;
        Likelihood cur_lik = target_->value();
        // before doing anything, calculate probability of sampling current state for adjustment
        Likelihood current_state_prob = calculateSamplingProb();

        // update infection time
        double cur_infection_time = infection_duration_->value();
        infection_duration_->saveState(stateId);

        const double eps           = normal_dist_(*rng_) * variance_;
        const double unconstrained = std::log(cur_infection_time - lower_bound_) - std::log(upper_bound_ - cur_infection_time);
        const double exp_prop      = std::exp(eps + unconstrained);
        const double proposed_infection_time = std::clamp((upper_bound_ * exp_prop + lower_bound_) / (exp_prop + 1), lower_bound_, upper_bound_);

        auto adj = log(proposed_infection_time - lower_bound_) +
                   log(upper_bound_ - proposed_infection_time) -
                   log(cur_infection_time - lower_bound_) -
                   log(upper_bound_ - cur_infection_time);

        infection_duration_->setValue(proposed_infection_time);

        // Genotypes are proposed from the parent set at the proposed infection time, as calculateSamplingProb() below
        // assumes; the reverse move would propose them from the parent set at the current time.
        genetics::proposeGenotypesFromParentSet<MaxParentSetSize>(infection_, latent_parent_, ps_->value(), *rng_, fn_rate_, fp_rate_, stateId);

        Likelihood proposed_state_prob = calculateSamplingProb();

        const auto acceptanceRatio = target_->value() - cur_lik + current_state_prob - proposed_state_prob + adj;
        // if (acceptanceRatio >= std::numeric_limits<double>::infinity()) {
        //     fmt::print("Acceptance ratio is infinity\n");
        //     fmt::print("Current state prob: {}\n", current_state_prob);
        //     fmt::print("Proposed state prob: {}\n", proposed_state_prob);
        //     fmt::print("Current Likelihood: {}\n", target_->value());
        //     fmt::print("Prev Likelihood: {}\n", cur_lik);
        //     fmt::print("Adj: {}\n", adj);
        // }
        const auto logProbAccept   = log(uniform_dist_(*rng_));
        const bool accept          = !std::isinf(acceptanceRatio) && logProbAccept <= acceptanceRatio;

        if (accept) {
            acceptances_++;
            countProposal(true);

            infection_duration_->acceptState();
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->acceptState();
            }
        } else {
            rejections_++;
            countProposal(false);
            infection_duration_->restoreState(stateId);
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->restoreState(stateId);
            }
        }

        assert(!target_->isDirty());

        total_updates_++;
    }


    template<typename T, typename Engine, typename InfectionEventImpl, typename GeneticsImpl, typename ParentSetImpl, int MaxParentSetSize, int MaxCOI>
    Likelihood JointGeneticsTimeSampler<T, Engine, InfectionEventImpl, GeneticsImpl, ParentSetImpl, MaxParentSetSize, MaxCOI>::calculateSamplingProb() noexcept {
        // The probability of proposing the current genotypes (see ParentSetGenotypeProposal.h): a uniform mixture
        // over the parent-set options.
        return proposalProbability_.logProb(*infection_, *latent_parent_, ps_->value());
    }

}

#endif//TRANSMISSION_NETWORKS_APP_JOINTGENETICSTIMESAMPLER_H
