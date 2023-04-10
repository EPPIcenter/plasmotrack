//
// Created by mmurphy on 10/25/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
#define TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H

#include "core/computation/transformers/Tempered.h"

#include <boost/random.hpp>
#include <fmt/core.h>

#include <filesystem>
#include <omp.h>

namespace fs = std::filesystem;

namespace transmission_nets::core::samplers {
    using namespace core::computation;
    // Implements the replication exchange (parallel tempering) algorithm

    template<typename State, typename Model, template<typename...> typename Scheduler, typename ModelLogger, typename StateLogger, typename Engine = boost::random::mt19937>
    struct ReplicaExchange {
//        using TemperedTarget = Tempered<Model>;
        struct Chain {
            std::shared_ptr<Model> target;
            std::shared_ptr<Model> model;
            std::shared_ptr<State> state;
            std::shared_ptr<Scheduler<Model>> sampler;
            std::shared_ptr<ModelLogger> modelLogger;
            std::shared_ptr<StateLogger> stateLogger;
            std::shared_ptr<Engine> r;
        };

        template<typename... Args>
        ReplicaExchange(int numChains, int samplesPerStep, double gradient, std::shared_ptr<Engine> r, fs::path outputDir, bool hotload, bool null_model, unsigned int num_cores, Args... args) : r_(r), num_cores_(num_cores) {
            swap_acceptance_rates = std::vector<int>(numChains, 0);
            for (int ii = 0; ii < numChains; ++ii) {
                swap_indices.push_back(ii);
                auto chain_r     = std::make_shared<Engine>(seed_dist_(*r_));
                double temp = 1.0 / std::pow(gradient, ii);

                std::shared_ptr<State> state;
                if (hotload) {
                    state = std::make_shared<State>(args..., chain_r, outputDir, null_model);
                } else {
                    state = std::make_shared<State>(args..., chain_r, null_model);
                }

                auto model       = std::make_shared<Model>(state, temp);
                auto target      = model;
//                auto target      = std::make_shared<TemperedTarget>(model, temp);
                auto sampler     = std::make_shared<Scheduler<Model>>(state, target, chain_r, samplesPerStep);
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
#pragma omp parallel default(none) num_threads(num_cores_)
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
                auto& chain_a = chains[swap_indices[ii]];
                auto& chain_b = chains[swap_indices[ii + 1]];

                fmt::print("Current llik: ({0:.3f}) {1:.2f} -- {2:.2f} ({3})\n", chain_a.target->getTemperature(), chain_a.target->value() / chain_a.target->getTemperature(), chain_a.target->value(), swap_indices[ii]);
                fmt::print("Current llik: ({0:.3f}) {1:.2f} -- {2:.2f} ({3})\n", chain_b.target->getTemperature(), chain_b.target->value() / chain_b.target->getTemperature(), chain_b.target->value(), swap_indices[ii + 1]);

                Likelihood curr_llik_a = chain_a.target->value();
                double temp_a          = chain_a.target->getTemperature();
                Likelihood curr_llik_b = chain_b.target->value();
                double temp_b          = chain_b.target->getTemperature();

                chain_a.target->setTemperature(temp_b);
                chain_b.target->setTemperature(temp_a);

                Likelihood prop_llik_a = chain_a.target->value();
                Likelihood prop_llik_b = chain_b.target->value();

                long double acceptanceRatio = (prop_llik_a - curr_llik_a + prop_llik_b - curr_llik_b);
                const bool accept      = log(uniform_dist_(*r_)) < acceptanceRatio;

                if (accept) {
                    swap_acceptance_rates[swap_indices[ii]]++;
                    swap_acceptance_rates[swap_indices[ii + 1]]++;
                    std::swap(swap_indices[ii], swap_indices[ii + 1]);
                    std::cout << fmt::format("Accepted Swap: {} (Hot: {}) {}", ii, swap_indices[0], hotValue()) << std::endl;
                } else {
                    chain_a.target->setTemperature(temp_a);
                    chain_b.target->setTemperature(temp_b);
                }
            }
            num_swaps++;
        }

        void logModel() {
            chains[swap_indices[0]].modelLogger->log();
        }

        void logState() {
            chains[swap_indices[0]].stateLogger->log();
        }

        Likelihood hotValue() {
            return chains[swap_indices[0]].target->value();
        }

        void printModelLlik() {
            fmt::print("Current P/L: {0:.2f} / {1:.2f}", chains[swap_indices[0]].model->getPrior(), chains[swap_indices[0]].model->getLikelihood());
        }



        void finalize() {
            chains[swap_indices[0]].stateLogger->finalize();
            chains[swap_indices[0]].modelLogger->finalize();
        }


        std::vector<Chain> chains{};
        std::vector<int> swap_acceptance_rates{};
        std::vector<size_t> swap_indices{};
        int num_swaps   = 0;
        std::shared_ptr<boost::random::mt19937> r_;
        unsigned int num_cores_;
        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<unsigned int> seed_dist_{};
    };


    //    template<typename T, typename State, typename Scheduler>
    //    ReplicaExchange<T, State, Scheduler>::ReplicaExchange(TemperedTarget &target, State &state, Scheduler &sampler) {
    //        addChain(target, state, sampler);
    //    }

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
