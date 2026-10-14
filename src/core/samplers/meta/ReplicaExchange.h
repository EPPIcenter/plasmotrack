//
// Created by mmurphy on 10/25/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
#define TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H

#include "core/samplers/meta/EnsembleSnapshot.h"
#include "spline.h"

#include <boost/random.hpp>

#include <fmt/core.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <ranges>
#include <stdexcept>
#include <vector>

// OpenMP support - conditionally include if available
#ifdef _OPENMP
#include <omp.h>
#elif __has_include(<omp.h>)
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace transmission_nets::core::samplers {
    using namespace core::computation;
    // Implements the replication exchange (parallel tempering) algorithm

    template<typename State, typename Model, template<typename...> typename Scheduler, typename Engine = boost::random::mt19937>
    struct ReplicaExchange {
        // The ladder adapts on the burnin() calls, numbered from 0, above LADDER_ADAPTS_AFTER that are multiples of
        // LADDER_ADAPTS_EVERY: 15, 20, ...
        static constexpr int LADDER_ADAPTS_AFTER = 10;
        static constexpr int LADDER_ADAPTS_EVERY = 5;

        struct Chain {
            std::shared_ptr<Model> target;
            std::shared_ptr<Model> model;
            std::shared_ptr<State> state;
            std::shared_ptr<Scheduler<Model>> sampler;
            std::shared_ptr<Engine> r;
        };

        // Every chain starts from a fresh State built from `args`.
        template<typename... Args>
        ReplicaExchange(const int numChains, const double gradient, std::shared_ptr<Engine> r, const bool null_model, const unsigned int num_cores, Args... args)
            : ReplicaExchange(numChains, gradient, std::move(r), num_cores, [=](size_t, std::shared_ptr<Engine> chain_r) {
                  return std::make_shared<State>(args..., chain_r, null_model);
              }) {}

        // Builds each replica's State with `makeState(replica, chainRng)` before its Model, so a State filled in by
        // `makeState` (as a checkpoint restore does) has every Model cache computed from it.
        using MakeState = std::function<std::shared_ptr<State>(size_t replica, std::shared_ptr<Engine> chainRng)>;

        ReplicaExchange(const int numChains, const double gradient, std::shared_ptr<Engine> r, const unsigned int num_cores, const MakeState& makeState) : r_(r), num_cores_(num_cores) {
            swap_counts.resize(numChains - 1);
            swap_barriers.resize(numChains - 1, 0.0);
            const double temp_step = (1.0 - gradient) / static_cast<double>(numChains);
            for (int ii = 0; ii < numChains; ++ii) {
                swap_indices.push_back(ii);
                auto chain_r = std::make_shared<Engine>(seed_dist_(*r_));
                double temp = 1.0 - temp_step * static_cast<double>(ii);
                temp_gradient.push_back(temp);

                std::shared_ptr<State> state = makeState(ii, chain_r);

                auto model = std::make_shared<Model>(state, temp);
                auto target = model;
                auto sampler = std::make_shared<Scheduler<Model>>(state, target, chain_r);

                Chain chain{target, model, state, sampler, chain_r};
                chains.push_back(chain);
            }
            for (const auto& chain : chains) {
                auto& replicaSamplers = samplers_.emplace_back();
                for (const auto& id : samplerIds()) {
                    replicaSamplers.push_back(&chain.sampler->scheduler_.sampler(id));
                }
            }
        }

        // Samplers adapt their proposals during burn-in only; adapting while sampling would change the chain's
        // stationary distribution.
        void sample() {
#pragma omp parallel default(none) num_threads(num_cores_)
            {
#pragma omp for
                for (size_t ii = 0; ii < chains.size(); ++ii) {
                    chains[ii].sampler->step(false);
                }
            }
            if (chains.size() > 1) {
                swapChains(false);
            }
        }

        // Returns the ladder adaptation this iteration made, if any: at burn-in iteration 15 (counting from 0), then every
        // 5. With two chains there is no interior temperature slot to move, so the ladder never adapts.
        std::optional<TemperatureAdaptation> burnin() {
#pragma omp parallel default(none) num_threads(num_cores_)
            {
#pragma omp for
                for (size_t ii = 0; ii < chains.size(); ++ii) {
                    chains[ii].sampler->step(true);
                }

            }
            if (chains.size() > 1) {
                swapChains(true);
            }
            const int iteration = burn_in_iterations++;
            if (chains.size() > 2 and iteration > LADDER_ADAPTS_AFTER and iteration % LADDER_ADAPTS_EVERY == 0) {
                return adaptTemperatures();
            }
            return std::nullopt;
        }

        void swapChains(const bool burnin) {
            for (size_t ii = even_swap; ii < chains.size() - 1; ii += 2) {
                auto& chain_a = chains[swap_indices[ii]];
                auto& chain_b = chains[swap_indices[ii + 1]];

                const Likelihood V_a = -chain_a.target->getLikelihood();
                const double temp_a = chain_a.target->getTemperature();
                const Likelihood V_b = -chain_b.target->getLikelihood();
                const double temp_b = chain_b.target->getTemperature();

                const double acceptance_ratio = (temp_b - temp_a) * (V_b - V_a);
                const double acceptance_rate = std::min(1.0, std::exp(acceptance_ratio));

                if (burnin) {
                    swap_barriers[ii] += 1.0 - acceptance_rate;
                }

                const double u = std::log(uniform_dist_(*r_));

                swap_counts[ii].attempted++;
                if ((acceptance_ratio > 0 || u < acceptance_ratio) and !std::isnan(acceptance_ratio)) {
                    std::swap(swap_indices[ii], swap_indices[ii + 1]);
                    chain_a.target->setTemperature(temp_b);
                    chain_b.target->setTemperature(temp_a);
                    swap_counts[ii].accepted++;
                }
            }

            if (burnin) {
                num_swaps++;
            }
            even_swap = !even_swap;
        }

        // The replica currently at the target posterior. Changes when chains swap.
        const Chain& hotChain() const {
            return chains[swap_indices[0]];
        }

        // What occupies each temperature slot right now.
        std::vector<TemperatureSlot> slots() const {
            std::vector<TemperatureSlot> view;
            view.reserve(chains.size());
            for (const size_t replica : swap_indices) {
                const auto& model = chains[replica].model;
                view.push_back({replica, model->getTemperature(), model->getPrior(), model->getLikelihood()});
            }
            return view;
        }

        // One entry per neighbouring pair of temperature slots; entry `pair` is slots `pair` and `pair + 1`.
        const std::vector<SwapCounts>& swapCounts() const {
            return swap_counts;
        }

        // The `<type>/<target>` ids of every replica's samplers, which all replicas register alike, in registration order.
        const std::vector<std::string>& samplerIds() const {
            return chains.front().sampler->scheduler_.samplerIds();
        }

        // What the Run recorder is handed after each Sampling iteration.
        EnsembleSnapshot<State, Model> snapshot() const {
            const auto& hot = hotChain();
            EnsembleSnapshot<State, Model> view{.hot = {*hot.state, *hot.model}, .slots = slots(), .swaps = swap_counts};
            for (const auto& replicaSamplers : samplers_) {
                auto& counts = view.acceptance.emplace_back();
                for (const auto* sampler : replicaSamplers) {
                    counts.push_back(sampler->counts());
                }
            }
            return view;
        }

        // Where adaptation stands: the Hot chain's samplers, by sampler id, and the ladder. Burn-in end fixes it.
        FinalAdaptation finalAdaptation() const {
            FinalAdaptation adaptation;
            for (const auto* sampler : samplers_[swap_indices[0]]) {
                adaptation.samplers.push_back(sampler->adaptationState());
            }
            for (const auto& slot : slots()) {
                adaptation.temperatures.push_back(slot.temperature);
            }
            return adaptation;
        }

        // Everything replica exchange carries between Sampling iterations, apart from each chain's own State and RNG.
        struct ExchangeState {
            std::vector<size_t> replicas;         // replica occupying each temperature slot
            std::vector<double> temperatures;     // ladder by temperature slot
            std::vector<double> swapBarriers;     // summed swap rejection per pair since the last adaptation
            int barrierIterations = 0;            // burn-in iterations those sums cover
            int burnInIterations  = 0;            // burn-in iterations run, which decide when the ladder next adapts
            bool secondPairSetNext = false;       // the next iteration tries pairs (1,2),(3,4),... instead of (0,1),...
            std::vector<SwapCounts> swapCounts;

            bool operator==(const ExchangeState&) const = default;
        };

        ExchangeState exchangeState() const {
            return {swap_indices, temp_gradient, swap_barriers, num_swaps, burn_in_iterations, even_swap, swap_counts};
        }

        // Puts each replica in its saved temperature slot at the saved temperature. Throws std::invalid_argument for a
        // state saved from an ensemble with a different number of chains.
        void restoreExchangeState(const ExchangeState& saved) {
            const size_t numSlots = chains.size();
            const size_t numPairs = numSlots - 1;
            if (saved.replicas.size() != numSlots or saved.temperatures.size() != numSlots or
                saved.swapBarriers.size() != numPairs or saved.swapCounts.size() != numPairs) {
                throw std::invalid_argument(fmt::format("exchange state does not fit an ensemble of {} chains", numSlots));
            }
            std::vector<bool> seen(numSlots, false);
            for (const size_t replica : saved.replicas) {
                if (replica >= numSlots or seen[replica]) {
                    throw std::invalid_argument("exchange state does not place each replica in exactly one temperature slot");
                }
                seen[replica] = true;
            }

            for (size_t slot = 0; slot < numSlots; ++slot) {
                chains[saved.replicas[slot]].target->setTemperature(saved.temperatures[slot]);
            }
            swap_indices  = saved.replicas;
            temp_gradient = saved.temperatures;
            swap_barriers = saved.swapBarriers;
            num_swaps     = saved.barrierIterations;
            burn_in_iterations = saved.burnInIterations;
            even_swap     = saved.secondPairSetNext;
            swap_counts   = saved.swapCounts;
        }

        Likelihood hotValue() {
            return chains[swap_indices[0]].target->value();
        }

        Likelihood hotValueThreadSafe() {
            return chains[swap_indices[0]].target->valueThreadSafe();
        }


        std::vector<Chain> chains{};

    private:
        // Respaces the Tempered chains' temperatures so that swap rejection is even across neighbouring slots. The Hot
        // chain stays at 1 and the lowest temperature stays put. A ladder that comes out NaN is not applied.
        TemperatureAdaptation adaptTemperatures() {
            const auto attempts = pairAttemptsSinceAdaptation();
            std::vector<double> barriers(swap_barriers.size(), 0.0);
            for (size_t ii = 0; ii < barriers.size(); ii++) {
                if (attempts[ii] > 0) {
                    barriers[ii] = swap_barriers[ii] / attempts[ii];
                }
            }
            // A pair with no attempts says nothing about its barrier; keep the ladder and the sums until it has some.
            if (std::ranges::find(attempts, 0) != attempts.end()) {
                return {false, temp_gradient, barriers};
            }

            // swap rate starts at t = 1 so need to reverse
            const std::vector<double> reversed_barriers(barriers.rbegin(), barriers.rend());

            std::vector<double> cumulative_swap_rate(chains.size(), 0.0);

            cumulative_swap_rate.front() = 0.0;
            for (size_t ii = 1; ii < cumulative_swap_rate.size(); ii++) {
                cumulative_swap_rate[ii] = cumulative_swap_rate[ii - 1] + reversed_barriers[ii - 1];
            }

            const std::vector<double> reversed_gradient{temp_gradient.rbegin(), temp_gradient.rend()};

            // monotonic cubic spline
            const tk::spline spl(reversed_gradient, cumulative_swap_rate, tk::spline::cspline, true);

            // target swap rates
            std::vector<double> cumulative_swap_grid(cumulative_swap_rate.size());
            cumulative_swap_grid.front() = cumulative_swap_rate.front();
            cumulative_swap_grid.back() = cumulative_swap_rate.back();
            const double step = (cumulative_swap_grid.back() - cumulative_swap_grid.front()) / (static_cast<double>(cumulative_swap_grid.size()) - 1.0);

            for (size_t ii = 1; ii < cumulative_swap_grid.size() - 1; ii++) {
                cumulative_swap_grid[ii] = cumulative_swap_grid[0] + static_cast<double>(ii) * step;
            }

            std::vector<double> new_temp_gradient(temp_gradient.size());
            new_temp_gradient.front() = temp_gradient.back();
            new_temp_gradient.back() = 1.0;
            for (size_t ii = 1; ii < temp_gradient.size() - 1; ii++) {
                new_temp_gradient[ii] = spl.solve_one(cumulative_swap_grid[ii]);
            }

            std::ranges::reverse(new_temp_gradient);

            if (std::ranges::any_of(new_temp_gradient, [](const double x) { return std::isnan(x); })) {
                return {false, temp_gradient, barriers};
            }

            for (size_t ii = 1; ii < temp_gradient.size() - 1; ii++) {
                chains[swap_indices[ii]].target->setTemperature(new_temp_gradient[ii]);
            }

            temp_gradient = new_temp_gradient;
            num_swaps = 0;
            std::ranges::fill(swap_barriers, 0.0);
            return {true, temp_gradient, barriers};
        }

        // Swap attempts per pair over the burn-in iterations since the last adaptation. Iterations alternate between
        // the pairs starting at slot 0 and those starting at slot 1, and the last one tried the set `even_swap` does
        // not name, so that set had the odd attempt out. Derived rather than counted, so a checkpoint needn't hold it.
        std::vector<int> pairAttemptsSinceAdaptation() const {
            const size_t lastSet = even_swap ? 0 : 1;
            std::vector<int> attempts(swap_barriers.size());
            for (size_t ii = 0; ii < attempts.size(); ii++) {
                attempts[ii] = ii % 2 == lastSet ? (num_swaps + 1) / 2 : num_swaps / 2;
            }
            return attempts;
        }

        // Read and restored through slots(), swapCounts() and exchangeState().
        std::vector<SwapCounts> swap_counts{};
        std::vector<std::vector<const AbstractSampler*>> samplers_{};// by replica, then by sampler id
        std::vector<size_t> swap_indices{};

        std::vector<double> swap_barriers{};
        std::vector<double> temp_gradient{};

        int num_swaps = 0;
        int burn_in_iterations = 0;// burnin() calls so far, which set when the ladder adapts
        bool even_swap = false;

        std::shared_ptr<boost::random::mt19937> r_;
        unsigned int num_cores_;
        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<unsigned int> seed_dist_{};
    };


}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
