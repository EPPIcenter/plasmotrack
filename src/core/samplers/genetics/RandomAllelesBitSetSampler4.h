//
// Created by Maxwell Murphy on 1/9/24.
//
#ifndef RANDOMALLELEBITSETSAMPLER4_H
#define RANDOMALLELEBITSETSAMPLER4_H

#include "core/containers/AllowedRelationships.h"
#include "core/datatypes/Data.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"

#include <boost/random.hpp>

#include <ranges>

#include <memory>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    class RandomAllelesBitSetSampler4 : public AbstractSampler {
    public:
        RandomAllelesBitSetSampler4(
                std::shared_ptr<InfectionEventImpl> infection,
                std::shared_ptr<InfectionEventImpl> latentParent,
                std::shared_ptr<LocusImpl> locus,
                std::shared_ptr<containers::AllowedRelationships<InfectionEventImpl>> allowedRelationships,
                std::shared_ptr<T> target,
                std::shared_ptr<Engine> rng,
                unsigned int max_coi) noexcept;

        void update() noexcept override;
        [[nodiscard]] unsigned int acceptances() noexcept;
        [[nodiscard]] unsigned int rejections() noexcept;
        [[nodiscard]] double acceptanceRate() noexcept;

    private:
        std::shared_ptr<InfectionEventImpl> infection_;
        std::shared_ptr<InfectionEventImpl> latentParent_;
        std::shared_ptr<LocusImpl> locus_;
        std::shared_ptr<containers::AllowedRelationships<InfectionEventImpl>> allowedRelationships_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int max_coi_;

        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> categorical_dist_{};

        unsigned int acceptances_ = 0;
        unsigned int rejections_ = 0;
        unsigned int total_updates_ = 0;

        // todo: make this tunable
        const double fp_rate = 0.05;
        const double fn_rate = 0.20;
        const double tx_from_parent = 0.90;
        const double tx_from_other = 0.20;
    };


    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, AlleleImpl, LocusImpl>::RandomAllelesBitSetSampler4(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<InfectionEventImpl> latentParent,
            std::shared_ptr<LocusImpl> locus,
            std::shared_ptr<containers::AllowedRelationships<InfectionEventImpl>> allowedRelationships,
            std::shared_ptr<T> target,
            std::shared_ptr<Engine> rng,
            unsigned int max_coi) noexcept
        : infection_(std::move(infection)),
          latentParent_(std::move(latentParent)),
          locus_(std::move(locus)),
          allowedRelationships_(allowedRelationships),
          target_(std::move(target)),
          rng_(std::move(rng)),
          max_coi_(max_coi){};

    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    void RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, AlleleImpl, LocusImpl>::update() noexcept {
        using Parameter = parameters::Parameter<AlleleImpl>;
        using Data = datatypes::Data<AlleleImpl>;
        namespace views = std::ranges::views;

        SAMPLER_STATE_ID stateId = RandomAlleleBitSet4ID;
        const Likelihood curLik = target_->value();

        std::vector<std::shared_ptr<Parameter>> updatedParameters{};

        // randomly select an allele state to be propagated through the transmission tree
        categorical_dist_.param(boost::random::uniform_int_distribution<>::param_type(0, locus_->totalAlleles() - 1));
        const int allele_idx = categorical_dist_(*rng_);

        // Log proposal probabilities of the proposed path given the current one, and of the current path given the
        // proposed one. Each infection's proposal depends on the allele in the infection before it on the path
        // (for the first, the latent parent's fixed anchor genotype), which the forward move has already changed and
        // the reverse move has not.
        double forward_log_prob = 0.0;
        double reverse_log_prob = 0.0;
        std::shared_ptr<InfectionEventImpl> curr_inf = infection_;
        bool prev_curr_state = latentParent_->latentGenotype(locus_)->value().allele(allele_idx) == 1;
        bool prev_proposed_state = prev_curr_state;
        int distance = 4;
        std::vector<int> path{};
        std::vector<int> curr_path{};
        std::vector<int> obs_path{};
        bool path_differs = false;
        while (distance > 0) {

            // collect potential children where infection time is greater than current infection
            std::vector<std::shared_ptr<InfectionEventImpl>> potential_children{};
            for (const auto& inf : allowedRelationships_->allowedChildren(curr_inf)) {
                if (inf->infectionTime() > curr_inf->infectionTime()) {
                    potential_children.push_back(inf);
                }
            }

            std::shared_ptr<Parameter> curr_param = curr_inf->latentGenotype(locus_);
            auto curr_val = curr_param->value();
            std::shared_ptr<Data> curr_data = nullptr;
            // if the locus is observed, then set curr data
            if (curr_inf->observedGenotype().contains(locus_)) {
                curr_data = curr_inf->observedGenotype(locus_);
            }

            // randomly select a child infection
            // (only when there is one: size() - 1 on an empty list wraps to an invalid range)
            int child_idx = 0;
            std::shared_ptr<InfectionEventImpl> child;
            if (!potential_children.empty()) {
                categorical_dist_.param(boost::random::uniform_int_distribution<>::param_type(0, static_cast<int>(potential_children.size()) - 1));
                child_idx = categorical_dist_(*rng_);
                child     = potential_children.at(child_idx);
            }

            // check if observed in curr infection
            const bool has_data = curr_data != nullptr;
            const bool in_curr = has_data && curr_inf->observedGenotype(locus_)->value().allele(allele_idx) == 1;

            // Check if observed in child infection
            const bool child_has_data = child != nullptr && child->observedGenotype().contains(locus_);
            const bool in_child = child_has_data && child->observedGenotype(locus_)->value().allele(allele_idx) == 1;

            const bool curr_state = curr_param->value().allele(allele_idx) == 1;

            const double obs_data_component = has_data ? (in_curr ? (1 - fp_rate) : (fn_rate)) : 1.0;
            const double child_data_component =
                    child_has_data ? in_child ? tx_from_parent * (1 - fn_rate) + (1 - tx_from_parent) * fp_rate : tx_from_parent * fn_rate + (1 - tx_from_parent) * (1 - fp_rate)
                                   : 1.0;

            // Probability of proposing the allele present, given whether it is present in the previous infection.
            auto prob_present = [&](const bool in_parent) {
                const double parent_prob_component = in_parent ? (tx_from_parent) : (tx_from_other);
                return parent_prob_component * obs_data_component * child_data_component;
            };
            const double forward_prob_present = prob_present(prev_proposed_state);
            const double reverse_prob_present = prob_present(prev_curr_state);

            const bool proposed_state = uniform_dist_(*rng_) < forward_prob_present;

            forward_log_prob += std::log(proposed_state ? forward_prob_present : 1 - forward_prob_present);
            reverse_log_prob += std::log(curr_state ? reverse_prob_present : 1 - reverse_prob_present);

            path.push_back(proposed_state);
            curr_path.push_back(curr_state);
            obs_path.push_back(has_data ? in_curr ? 1 : 0 : -1);

            if (proposed_state != curr_state) {
                path_differs = true;
                curr_param->saveState(stateId);
                updatedParameters.push_back(curr_param);

                curr_val.set(allele_idx, proposed_state);
                curr_param->setValue(curr_val);
            }

            if (potential_children.empty()) {
                break;
            }

            prev_curr_state     = curr_state;
            prev_proposed_state = proposed_state;
            curr_inf = potential_children[child_idx];
            distance--;
        }

        if (!path_differs) {
            acceptances_++;
            countProposal(true);
            return;
        }

        const Likelihood adj = reverse_log_prob - forward_log_prob;
        const Likelihood acceptanceRatio = target_->value() - curLik + adj;
        const Likelihood logProbAccept = std::log(uniform_dist_(*rng_));
        const bool accept = logProbAccept <= acceptanceRatio;

        // if (debug_ and accept) {
        //     fmt::print("Acceptance ratio: {}\n", acceptanceRatio);
        //     fmt::print("Log probability of acceptance: {}\n", logProbAccept);
        //     fmt::print("Current likelihood: {}\n", curLik);
        //     fmt::print("Forward Proposal Log Probability: {}\n", forward_log_prob);
        //     fmt::print("Reverse Proposal Log Probability: {}\n", reverse_log_prob);
        //     fmt::print("Proposed likelihood: {}\n", target_->value());
        //     fmt::print("Adjustment: {}\n", adj);
        //     fmt::print("Proposed Path: {}\n", io::serialize(path));
        //     fmt::print("Current Path: {}\n", io::serialize(curr_path));
        //     fmt::print("Observed Path: {}\n", io::serialize(obs_path));
        //     fmt::print("----------------------------------------------------\n");
        // }

        if (accept) {
            acceptances_++;
            countProposal(true);
            for (auto& param : updatedParameters) {
                param->acceptState();
            }
        } else {
            rejections_++;
            countProposal(false);
            for (auto& param : updatedParameters) {
                param->restoreState(stateId);
            }
        }

        total_updates_++;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    unsigned int RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, AlleleImpl, LocusImpl>::acceptances() noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    unsigned int RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, AlleleImpl, LocusImpl>::rejections() noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename InfectionEventImpl, typename AlleleImpl, typename LocusImpl>
    double RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, AlleleImpl, LocusImpl>::acceptanceRate() noexcept {
        return static_cast<double>(acceptances_) / static_cast<double>(total_updates_);
    }


}// namespace transmission_nets::core::samplers::genetics

#endif//RANDOMALLELEBITSETSAMPLER4_H
