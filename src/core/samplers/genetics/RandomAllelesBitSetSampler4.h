//
// Created by Maxwell Murphy on 1/9/24.
//
#ifndef RANDOMALLELEBITSETSAMPLER4_H
#define RANDOMALLELEBITSETSAMPLER4_H

#include "core/datatypes/Genotype.h"
#include "core/containers/AllowedRelationships.h"
#include "core/datatypes/Data.h"
#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/genetics/ObservationRate.h"

#include <boost/random.hpp>

#include <algorithm>
#include <ranges>

#include <memory>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    class RandomAllelesBitSetSampler4 : public AbstractSampler {
    public:
        // `expectedFalsePositives` and `expectedFalseNegatives`, when set, are that infection's sampled counts. The
        // per-allele rates are those counts divided by this locus's allele count, read at each update. Without them
        // the proposal uses the fixed rates below.
        RandomAllelesBitSetSampler4(
                std::shared_ptr<InfectionEventImpl> infection,
                std::shared_ptr<InfectionEventImpl> latentParent,
                std::shared_ptr<LocusImpl> locus,
                std::shared_ptr<containers::AllowedRelationships<InfectionEventImpl>> allowedRelationships,
                std::shared_ptr<T> target,
                std::shared_ptr<Engine> rng,
                unsigned int max_coi,
                std::shared_ptr<parameters::Parameter<double>> expectedFalsePositives = nullptr,
                std::shared_ptr<parameters::Parameter<double>> expectedFalseNegatives = nullptr) noexcept;

        void update() noexcept override;
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

        // Used when the infection's sampled error counts are not supplied.
        double fp_rate = 0.05;
        double fn_rate = 0.20;
        const double tx_from_parent = 0.90;
        const double tx_from_other = 0.20;
        std::shared_ptr<parameters::Parameter<double>> expectedFalsePositives_;
        std::shared_ptr<parameters::Parameter<double>> expectedFalseNegatives_;
    };


    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, LocusImpl>::RandomAllelesBitSetSampler4(
            std::shared_ptr<InfectionEventImpl> infection,
            std::shared_ptr<InfectionEventImpl> latentParent,
            std::shared_ptr<LocusImpl> locus,
            std::shared_ptr<containers::AllowedRelationships<InfectionEventImpl>> allowedRelationships,
            std::shared_ptr<T> target,
                std::shared_ptr<Engine> rng,
                unsigned int max_coi,
                std::shared_ptr<parameters::Parameter<double>> expectedFalsePositives,
                std::shared_ptr<parameters::Parameter<double>> expectedFalseNegatives) noexcept
        : infection_(std::move(infection)),
          latentParent_(std::move(latentParent)),
          locus_(std::move(locus)),
          allowedRelationships_(allowedRelationships),
          target_(std::move(target)),
          rng_(std::move(rng)),
          max_coi_(max_coi),
          expectedFalsePositives_(std::move(expectedFalsePositives)),
          expectedFalseNegatives_(std::move(expectedFalseNegatives)){};

    template<typename T, typename Engine, typename InfectionEventImpl, typename LocusImpl>
    void RandomAllelesBitSetSampler4<T, Engine, InfectionEventImpl, LocusImpl>::update() noexcept {
        using Parameter = parameters::Parameter<core::datatypes::Genotype>;
        using Genotype = core::datatypes::Genotype;
        using Data = datatypes::Data<core::datatypes::Genotype>;
        namespace views = std::ranges::views;

        // The observation likelihood uses expected count / alleles at this locus. The same rates in both directions,
        // since this move does not change them.
        if (expectedFalsePositives_ and expectedFalseNegatives_ and locus_->totalAlleles() > 0) {
            const double alleles = static_cast<double>(locus_->totalAlleles());
            fp_rate = perAlleleRate(expectedFalsePositives_->value(), alleles);
            fn_rate = perAlleleRate(expectedFalseNegatives_->value(), alleles);
        }
        // randomly select an allele state to be propagated through the transmission tree
        categorical_dist_.param(boost::random::uniform_int_distribution<>::param_type(0, locus_->totalAlleles() - 1));
        const int allele_idx = categorical_dist_(*rng_);

        // Log proposal probabilities of the proposed path given the current one, and of the current path given the
        // proposed one. Each infection's proposal depends on the allele in the infection before it on the path
        // (for the first, the latent parent's fixed anchor genotype), which the forward move has already changed and
        // the reverse move has not.
        double forward_log_prob = 0.0;
        double reverse_log_prob = 0.0;
        // Applied only if some allele changes. A path that proposes the current alleles is not a proposal.
        std::vector<std::pair<std::shared_ptr<Parameter>, Genotype>> changes;
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
                curr_val.set(allele_idx, proposed_state);
                changes.emplace_back(curr_param, curr_val);
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
            return;
        }
        auto proposal = this->propose(*target_, *rng_);
        for (const auto& [param, value] : changes) {
            proposal.change(*param, value);
        }
        proposal.metropolis(reverse_log_prob - forward_log_prob);
    }


}// namespace transmission_nets::core::samplers::genetics

#endif//RANDOMALLELEBITSETSAMPLER4_H
