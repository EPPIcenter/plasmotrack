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
#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/generators/RandomSequence.h"
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
                double fn_rate = 0.05,
                double fp_rate = 0.05
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
        boost::random::uniform_int_distribution<> parent_set_sampling_dist{};

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
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
        core::utils::generators::CombinationIndicesGenerator ps_idx_gen;

        int total_possible_parent_sets = 0;
        std::vector<int> cumulative_possible_parent_set_sizes{};

        auto tmp_ps = ps_->value();

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

        // Calculate the total number of possible parent sets across possible parent set sizes
        for (size_t i = 1; i <= MaxParentSetSize and i <= tmp_ps.size(); ++i) {
            total_possible_parent_sets += boost::math::binomial_coefficient<double>(tmp_ps.size(), i);
            cumulative_possible_parent_set_sizes.push_back(total_possible_parent_sets);
        }

        bool include_latent_parent = false;
        std::vector<unsigned int> parent_set_idxs{};

        if (total_possible_parent_sets > 0) {
            parent_set_sampling_dist.param(boost::random::uniform_int_distribution<>::param_type(0, total_possible_parent_sets - 1));
            int parent_set_idx = parent_set_sampling_dist(*rng_);
            include_latent_parent = uniform_dist_(*rng_) < 0.5;

            // after sampling the index, advance the generator to the correct parent set
            for (size_t i = 0; i < cumulative_possible_parent_set_sizes.size(); ++i) {
                if (parent_set_idx < cumulative_possible_parent_set_sizes[i]) {
                    ps_idx_gen.reset(tmp_ps.size(), i + 1);
                    ps_idx_gen.advance(cumulative_possible_parent_set_sizes[i] - parent_set_idx);
                    break;
                }
            }
            parent_set_idxs = ps_idx_gen.curr;

        } else {
            // if there are no possible parent sets, then the latent parent must be included
            include_latent_parent = true;
        }

        // Sample a new proposed genetic state with at least one coming from each parent
        for (const auto& locus : infection_->loci()) {
            infection_->latentGenotype(locus)->saveState(stateId);
            auto proposal = infection_->latentGenotype(locus)->value();
            proposal.reset();

            // accumulate the shared alleles across all parents
            GeneticsImpl all_shared = latent_parent_->latentGenotype(locus)->value();
            if (!parent_set_idxs.empty()) {
                if (!include_latent_parent) {
                    all_shared = tmp_ps.begin()[parent_set_idxs[0]]->latentGenotype(locus)->value();
                    for (size_t i = 1; i < parent_set_idxs.size(); ++i) {
                        all_shared = GeneticsImpl::any(all_shared, tmp_ps.begin()[parent_set_idxs[i]]->latentGenotype(locus)->value());
                    }
                } else {
                    for (unsigned int& parent_set_idx : parent_set_idxs) {
                        all_shared = GeneticsImpl::any(all_shared, tmp_ps.begin()[parent_set_idx]->latentGenotype(locus)->value());
                    }
                }
            }

            const auto& rand_seq = core::utils::generators::randomSequence(0, proposal.totalAlleles(), rng_);

            // flag to indicate if we have set at least 1 allele
            boost::container::flat_map<int, bool> one_set_flags{};
            for (const auto& idx : parent_set_idxs) {
                one_set_flags[int(idx)] = false;
            }

            if (include_latent_parent) {
                one_set_flags[-1] = false;
            }


            if (infection_->observedGenotype().contains(locus)) {
                const auto& child_observed_genotype = infection_->observedGenotype(locus)->value();

                // for each allele, check if it's present in the parent set, then sample a new allele
                // conditional on the observed data. Make sure at least one allele is present from each parent
                for (const auto& i : rand_seq) {

                    std::vector<int> possible_parent_idxs{};
                    for (const auto& parent_set_idx : parent_set_idxs) {
                        if (tmp_ps.begin()[parent_set_idx]->latentGenotype(locus)->value().allele(i) and not one_set_flags[int(parent_set_idx)]) {
                            possible_parent_idxs.push_back(int(parent_set_idx));
                        }
                    }

                    if (include_latent_parent) {
                        if (latent_parent_->latentGenotype(locus)->value().allele(i) and not one_set_flags[-1]) {
                            possible_parent_idxs.push_back(-1);
                        }
                    }

                    if (all_shared.allele(i)) {
                        const auto p = uniform_dist_(*rng_);
                        // the allele is observed, so set to 1 with probability 1 - fp_rate
                        if (child_observed_genotype.allele(i)) {
                            if (p < 1 - fp_rate_ or not possible_parent_idxs.empty()) {
                                proposal.set(i, true);
                                for (const auto parent_set_idx : possible_parent_idxs) {
                                    one_set_flags[parent_set_idx] = true;
                                }
                            } else {
                                proposal.set(i, false);
                            }

                        } else {
                            // the allele is not observed, so set to 0 with probability 1 - fn_rate
                            if (p < 1 - fn_rate_ and possible_parent_idxs.empty()) {
                                proposal.set(i, false);
                            } else {
                                proposal.set(i, true);
                                for (const auto parent_set_idx : possible_parent_idxs) {
                                    one_set_flags[parent_set_idx] = true;
                                }
                            }
                        }
                    }
                }
                infection_->latentGenotype(locus)->setValue(proposal);
            } else {
                // if the locus is not observed, then sample a new allele for each allele in the parent set
                for (const auto& i : rand_seq) {

                    std::vector<int> possible_parent_idxs{};
                    for (const auto& parent_set_idx : parent_set_idxs) {
                        if (tmp_ps.begin()[parent_set_idx]->latentGenotype(locus)->value().allele(i) and not one_set_flags[int(parent_set_idx)]) {
                            possible_parent_idxs.push_back(int(parent_set_idx));
                        }
                    }

                    if (include_latent_parent) {
                        if (latent_parent_->latentGenotype(locus)->value().allele(i) and not one_set_flags[-1]) {
                            possible_parent_idxs.push_back(-1);
                        }
                    }

                    if (all_shared.allele(i)) {
                        const auto p = uniform_dist_(*rng_);
                        if (p < 0.5 or not possible_parent_idxs.empty()) {
                            proposal.set(i, true);
                            for (const auto& parent_set_idx : possible_parent_idxs) {
                                one_set_flags[parent_set_idx] = true;
                            }
                        } else {
                            proposal.set(i, false);
                        }
                    }
                }
                infection_->latentGenotype(locus)->setValue(proposal);
            }
        }

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

            infection_duration_->acceptState();
            for (const auto& locus : infection_->loci()) {
                infection_->latentGenotype(locus)->acceptState();
            }
        } else {
            rejections_++;
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
        // Calculate the probability of sampling the currently set alleles conditional
        // on possible parent sets and the observed genetics. An allele may be present
        // *only* if it is present in a parent set. If the allele is present in the observed
        // data, then it is present in the child with probability 1 - fp_rate, and fp_rate
        // otherwise. If the allele is not present in the observed data, then it is present
        // in the child with probability fn_rate, and 1 - fn_rate otherwise.
        std::vector<Likelihood> probs{};
        auto tmp_ps = ps_->value();
        const size_t total_parents = tmp_ps.size();
        size_t total_parent_sets = 0;
        Likelihood max_llik = -std::numeric_limits<Likelihood>::infinity();

        // todo: need to make sure we're accounting for some of the alleles being present w/ probability = 1 due to
        //  requiring every parent to propagate at least one allele
        // iterate over all possible parent sets
        for (size_t num_parents = 1; num_parents <= MaxParentSetSize and num_parents <= total_parents; ++num_parents) {
            core::utils::generators::CombinationIndicesGenerator psIdxGen(total_parents, num_parents);

            const auto& parent_set_idxs = psIdxGen.curr;
            // iterate over all possible parent sets of size num_parents
            while(!psIdxGen.completed) {

                // iterate over the case with and without the latent parent
                for (size_t ii = 0; ii < 2; ++ii) {
                    // under a particular parent set, how many alleles are tn, fn, fp, tp?
                    int total_tn = 0;
                    int total_fn = 0;
                    int total_fp = 0;
                    int total_tp = 0;
                    bool is_valid_parent_set = true;
                    for (const auto& locus : infection_->loci()) {
                        if (infection_->observedGenotype().contains(locus)) {
                            const auto& child_latent_genotype   = infection_->latentGenotype(locus)->value();
                            const auto& child_observed_genotype = infection_->observedGenotype(locus)->value();

                            // accumulate the shared alleles across all parents

                            GeneticsImpl all_shared = latent_parent_->latentGenotype(locus)->value();
                            if (!parent_set_idxs.empty()) {
                                if (ii != 0) {
                                    all_shared = tmp_ps.begin()[parent_set_idxs[0]]->latentGenotype(locus)->value();
                                    for (size_t i = 1; i < parent_set_idxs.size(); ++i) {
                                        all_shared = GeneticsImpl::any(all_shared, tmp_ps.begin()[parent_set_idxs[i]]->latentGenotype(locus)->value());
                                    }
                                } else {
                                    for (const unsigned int& parent_set_idx : parent_set_idxs) {
                                        all_shared = GeneticsImpl::any(all_shared, tmp_ps.begin()[parent_set_idx]->latentGenotype(locus)->value());
                                    }
                                }
                            }

                            int false_positives = GeneticsImpl::falsePositiveCount(all_shared, child_latent_genotype);

                            // the latent genotype is incompatible with the parent set if there are alleles in the child
                            // that are not in the parent set
                            if (false_positives > 0) {
                                probs.push_back(-std::numeric_limits<Likelihood>::infinity());

                                if (ii != 0) {
                                    psIdxGen.next();
                                    total_parent_sets++;
                                }
                                is_valid_parent_set = false;
                                break;
                            }

                            int illegal_fp_count = GeneticsImpl::falsePositiveCount(all_shared, child_observed_genotype);
                            int illegal_tn_count = GeneticsImpl::trueNegativeCount(all_shared, child_observed_genotype);
                            int fn_count         = GeneticsImpl::falseNegativeCount(child_latent_genotype, child_observed_genotype);
                            int tp_count         = GeneticsImpl::truePositiveCount(child_latent_genotype, child_observed_genotype);
                            int tn_count         = GeneticsImpl::trueNegativeCount(child_latent_genotype, child_observed_genotype) - illegal_tn_count;
                            int fp_count         = GeneticsImpl::falsePositiveCount(child_latent_genotype, child_observed_genotype) - illegal_fp_count;
                            assert((int) all_shared.totalPositiveCount() == tp_count + fn_count + tn_count + fp_count);

                            total_tn += tn_count;
                            total_fn += fn_count;
                            total_fp += fp_count;
                            total_tp += tp_count;
                        }
                    }
                    if (is_valid_parent_set) {
                        probs.push_back(total_tn * std::log(1 - fn_rate_) + total_fn * std::log(fn_rate_) + total_fp * std::log(fp_rate_) + total_tp * std::log(1 - fp_rate_));
                        max_llik = std::max(max_llik, probs.back());
                        if (ii != 0) {
                            psIdxGen.next();
                            total_parent_sets++;
                        }
                    }
                }
            }
        }


        // Include the latent parent as a single parent set
        int total_tn = 0;
        int total_fn = 0;
        int total_fp = 0;
        int total_tp = 0;
        bool is_valid_parent_set = true;
        for (const auto& locus : infection_->loci()) {
            if (infection_->observedGenotype().contains(locus)) {
                const auto& child_latent_genotype   = infection_->latentGenotype(locus)->value();
                const auto& child_observed_genotype = infection_->observedGenotype(locus)->value();

                GeneticsImpl all_shared = latent_parent_->latentGenotype(locus)->value();

                int false_positives = GeneticsImpl::falsePositiveCount(all_shared, child_latent_genotype);

                // the latent genotype is incompatible with the parent set if there are alleles in the child
                // that are not in the parent set
                if (false_positives > 0) {
                    probs.push_back(-std::numeric_limits<Likelihood>::infinity());
                    total_parent_sets++;
                    is_valid_parent_set = false;
                    break;
                }

                int illegal_fp_count = GeneticsImpl::falsePositiveCount(all_shared, child_observed_genotype);
                int illegal_tn_count = GeneticsImpl::trueNegativeCount(all_shared, child_observed_genotype);
                int fn_count         = GeneticsImpl::falseNegativeCount(child_latent_genotype, child_observed_genotype);
                int tp_count         = GeneticsImpl::truePositiveCount(child_latent_genotype, child_observed_genotype);
                int tn_count         = GeneticsImpl::trueNegativeCount(child_latent_genotype, child_observed_genotype) - illegal_tn_count;
                int fp_count         = GeneticsImpl::falsePositiveCount(child_latent_genotype, child_observed_genotype) - illegal_fp_count;
                assert((int) all_shared.totalPositiveCount() == tp_count + fn_count + tn_count + fp_count);

                total_tn += tn_count;
                total_fn += fn_count;
                total_fp += fp_count;
                total_tp += tp_count;

            }
        }
        if (is_valid_parent_set) {
            probs.push_back(total_tn * std::log(1 - fn_rate_) + total_fn * std::log(fn_rate_) + total_fp * std::log(fp_rate_) + total_tp * std::log(1 - fp_rate_));
            max_llik = std::max(max_llik, probs.back());
            total_parent_sets++;
        }


        return core::utils::logSumExpKnownMax(probs.begin(), probs.end(), max_llik) - std::log(total_parent_sets * 2 + 1);
    }

}

#endif//TRANSMISSION_NETWORKS_APP_JOINTGENETICSTIMESAMPLER_H
