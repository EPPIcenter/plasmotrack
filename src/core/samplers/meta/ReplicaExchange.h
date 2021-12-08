//
// Created by mmurphy on 10/25/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
#define TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H

#include "core/computation/transformers/Tempered.h"

#include <boost/random.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <omp.h>

namespace fs = std::filesystem;

namespace transmission_nets::core::samplers {
    using namespace core::computation;
    // Implements the replication exchange (parallel tempering) algorithm

    template<typename State, typename Model, template<typename...> typename Scheduler, typename ModelLogger, typename StateLogger, typename Engine = boost::random::mt19937>
    struct ReplicaExchange {
        using TemperedTarget = Tempered<Model>;
        struct Chain {
            std::shared_ptr<TemperedTarget> target;
            std::shared_ptr<Model> model;
            std::shared_ptr<State> state;
            std::shared_ptr<Scheduler<TemperedTarget>> sampler;
            std::shared_ptr<ModelLogger> modelLogger;
            std::shared_ptr<StateLogger> stateLogger;
            std::shared_ptr<Engine> r;
        };

        template<typename... Args>
        ReplicaExchange(int numChains, int samplesPerStep, double gradient, std::shared_ptr<Engine> r, fs::path outputDir, bool hotload, Args... args) : r_(r) {
            swap_acceptance_rates = std::vector<int>(numChains, 0);
            for (int ii = 0; ii < numChains; ++ii) {
                double temp = 1.0 / std::pow(gradient, ii);
                std::shared_ptr<State> state;
                if (hotload) {
                    state = std::make_shared<State>(args..., outputDir);
                } else {
                    state = std::make_shared<State>(args...);
                }
                auto model = std::make_shared<Model>(state);
                auto target = std::make_shared<TemperedTarget>(model, temp);
                auto chain_r = std::make_shared<Engine>(seed_dist_(*r_));
                auto sampler = std::make_shared<Scheduler<TemperedTarget>>(state, target, chain_r, samplesPerStep);
                auto modelLogger = std::make_shared<ModelLogger>(model, outputDir);
                std::shared_ptr<StateLogger> stateLogger;
                if (ii == (numChains - 1)) {
                    stateLogger = std::make_shared<StateLogger>(state, outputDir, true);
                } else {
                    stateLogger = std::make_shared<StateLogger>(state, outputDir, false);
                }
                Chain chain{target, model, state, sampler, modelLogger, stateLogger, chain_r};
                chains.push_back(chain);
            }
        }

        void sample() {

#pragma omp parallel default(none) num_threads(chains.size())
            {
#pragma omp for
                for (size_t ii = 0; ii < chains.size(); ++ii) {
                    chains[ii].sampler->step();
                }
            }
            swap_chains();
        }

        void swap_chains() {
            // Generate sequences of (0,2,4...) and (1,3,5...) every other swap
            for (size_t ii = 0 + num_swaps % 2; ii < chains.size() - 1; ii += 2) {
                fmt::print("Current llik: {0:.2f} -- {1:.2f} ({2})\n", chains[ii].target->value() / chains[ii].target->getTemperature(), chains[ii].target->value(), ii);
                fmt::print("Current llik: {0:.2f} -- {1:.2f} ({2})\n", chains[ii + 1].target->value() / chains[ii + 1].target->getTemperature(), chains[ii + 1].target->value(), ii + 1);
                Likelihood curr_llik_a = chains[ii].target->value();
                double temp_a = chains[ii].target->getTemperature();
                Likelihood curr_llik_b = chains[ii + 1].target->value();
                double temp_b = chains[ii + 1].target->getTemperature();

                chains[ii].target->setTemperature(temp_b);
                chains[ii + 1].target->setTemperature(temp_a);
                Likelihood prop_llik_a = chains[ii].target->value();
                Likelihood prop_llik_b = chains[ii + 1].target->value();

                double acceptanceRatio = (prop_llik_a - curr_llik_a + prop_llik_b - curr_llik_b);
                const bool accept = log(uniform_dist_(*r_)) < acceptanceRatio;

                if (accept) {
                    swap_acceptance_rates[ii]++;
                    swap_acceptance_rates[ii + 1]++;
                    if (ii == hot_idx_) {
                        hot_idx_++;
                    } else if ((ii + 1) == hot_idx_) {
                        hot_idx_--;
                    }
                    std::cout << fmt::format("Accepted Swap: {} (Hot: {}) {}", ii, hot_idx_, chains[hot_idx_].target->value()) << std::endl;
                } else {
                    chains[ii].target->setTemperature(temp_a);
                    chains[ii + 1].target->setTemperature(temp_b);
                }
            }
            num_swaps++;
        }

        void logModel() {
            chains[hot_idx_].modelLogger->log();
        }

        void logState() {
            chains[hot_idx_].stateLogger->log();
        }

        Likelihood hotValue() {
            return chains[hot_idx_].target->value();
        }

        void finalize() {
            chains[hot_idx_].stateLogger->finalize();
            chains[hot_idx_].modelLogger->finalize();
        }


        std::vector<Chain> chains{};
        std::vector<int> swap_acceptance_rates{};
        size_t hot_idx_ = 0;
        int num_swaps = 0;
        std::shared_ptr<boost::random::mt19937> r_;
        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<unsigned int> seed_dist_{};
    };


    //    template<typename T, typename State, typename Scheduler>
    //    ReplicaExchange<T, State, Scheduler>::ReplicaExchange(TemperedTarget &target, State &state, Scheduler &sampler) {
    //        addChain(target, state, sampler);
    //    }

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
