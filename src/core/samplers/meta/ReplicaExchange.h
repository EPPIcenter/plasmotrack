//
// Created by mmurphy on 10/25/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
#define TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H

#include "core/computation/transformers/Tempered.h"

#include "spline.h"

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
            swap_acceptance_rates.resize(numChains - 1, 0);
            swap_barriers.resize(numChains - 1, 0.0);
            double temp_step = (1.0 - gradient) / (double) (numChains);
            for (int ii = 0; ii < numChains; ++ii) {
                swap_indices.push_back(ii);
                auto chain_r = std::make_shared<Engine>(seed_dist_(*r_));
                double temp = 1.0 - temp_step * ii;
                temp_gradient.push_back(temp);

                std::shared_ptr<State> state;
                if (hotload) {
                    state = std::make_shared<State>(args..., chain_r, outputDir, null_model);
                } else {
                    state = std::make_shared<State>(args..., chain_r, null_model);
                }

                auto model  = std::make_shared<Model>(state, temp);
                auto target = model;
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
            if (chains.size() > 1) {
                swapChains(false, false);
            }
//            swap_chains(false, false);
        }

        void burnin() {
#pragma omp parallel default(none) num_threads(num_cores_)
            {
#pragma omp for
                for (size_t ii = 0; ii < chains.size(); ++ii) {
                    chains[ii].sampler->step();
                }
            }
            if (chains.size() > 1) {
                swapChains(true, false);
            }
//            swap_chains(true, false);
        }

//        void swap_chains(bool burnin, bool pre_adapt_temperature) {
//            // Generate sequences of (0,2,4...) and (1,3,5...) every other swap
//            for (size_t ii = 0 + num_swaps % 2; ii < chains.size() - 1; ii += 2) {
//                auto& chain_a = chains[swap_indices[ii]];
//                auto& chain_b = chains[swap_indices[ii + 1]];
//
//                fmt::print("Current llik: ({-1:.3f}) {1:.2f} -- {2:.2f} ({3})\n", chain_a.target->getTemperature(), chain_a.target->value() / chain_a.target->getTemperature(), chain_a.target->value(), swap_indices[ii]);
//                fmt::print("Current llik: ({-1:.3f}) {1:.2f} -- {2:.2f} ({3})\n", chain_b.target->getTemperature(), chain_b.target->value() / chain_b.target->getTemperature(), chain_b.target->value(), swap_indices[ii + 1]);
//
//                Likelihood curr_llik_a = chain_a.target->value();
//                double temp_a          = chain_a.target->getTemperature();
//                Likelihood curr_llik_b = chain_b.target->value();
//                double temp_b          = chain_b.target->getTemperature();
//
//                chain_a.target->setTemperature(temp_b);
//                chain_b.target->setTemperature(temp_a);
//
//                Likelihood prop_llik_a = chain_a.target->value();
//                Likelihood prop_llik_b = chain_b.target->value();
//
//                long double acceptanceRatio = (prop_llik_a - curr_llik_a + prop_llik_b - curr_llik_b);
//                const bool accept           = log(uniform_dist_(*r_)) < acceptanceRatio;
//
//                if (accept) {
//                    swap_acceptance_rates[swap_indices[ii]]++;
//                    swap_acceptance_rates[swap_indices[ii + 1]]++;
//                    std::swap(swap_indices[ii], swap_indices[ii + 1]);
//                    std::cout << fmt::format("Accepted Swap: {} (Hot: {}) {}", ii, swap_indices[0], hotValue()) << std::endl;
//                } else {
//                    chain_a.target->setTemperature(temp_a);
//                    chain_b.target->setTemperature(temp_b);
//                }
//            }
//            num_swaps++;
//        }

        void swapChains(bool burnin, bool pre_adapt_temperature) {
            for (size_t ii = (int) even_swap; ii < chains.size() - 1; ii += 2) {
                auto& chain_a = chains[swap_indices[ii]];
                auto& chain_b = chains[swap_indices[ii + 1]];
                Likelihood V_a = -chain_a.target->getLikelihood();
                double temp_a  = chain_a.target->getTemperature();

                Likelihood V_b = -chain_b.target->getLikelihood();
                double temp_b  = chain_b.target->getTemperature();

                long double acceptance_ratio = (temp_b - temp_a) * (V_b - V_a);

                double acceptance_rate = std::min(1.0, (double) std::exp(acceptance_ratio));

                if (burnin and !pre_adapt_temperature) {
                    swap_barriers[ii] += 1.0 - acceptance_rate;
                }

                double u = std::log(uniform_dist_(*r_));

                if ((acceptance_ratio > 0 || u < acceptance_ratio) and !std::isnan(acceptance_ratio)) {
                    std::swap(swap_indices[ii], swap_indices[ii + 1]);
                    chain_a.target->setTemperature(temp_b);
                    chain_b.target->setTemperature(temp_a);

                    if (!burnin) {
                        swap_acceptance_rates[ii]++;
                    }
                    fmt::print("Accepted Swap: {} (Hot: {}) {}\n", ii, swap_indices[0], hotValue());
                }
            }
            swap_store.push_back(swap_indices[0]);

            if (burnin and !pre_adapt_temperature) {
                num_swaps++;
            }
            even_swap = !even_swap;
        }

        void adaptTemp() {

            // swap rate starts at t = 1 so need to reverse
            const std::vector<double> reversed_swap_barriers(swap_barriers.rbegin(), swap_barriers.rend());

            std::vector<double> cumulative_swap_rate = std::vector<double>(chains.size(), 0.0);

            cumulative_swap_rate.front() = 0.0;
            for (size_t ii = 1; ii < cumulative_swap_rate.size(); ii++) {
                cumulative_swap_rate[ii] = cumulative_swap_rate[ii - 1] + reversed_swap_barriers[ii - 1] / ((double) num_swaps / 2);
            }

            const std::vector<double> reversed_gradient{temp_gradient.rbegin(), temp_gradient.rend()};

            // monotonic cubic spline
            const tk::spline spl(reversed_gradient, cumulative_swap_rate, tk::spline::cspline, true);

            // target swap rates
            std::vector<double> cumulative_swap_grid              = std::vector<double>(cumulative_swap_rate.size());
            cumulative_swap_grid.front()                          = cumulative_swap_rate.front();
            cumulative_swap_grid[cumulative_swap_grid.size() - 1] = cumulative_swap_rate.back();
            const double step                                     = (cumulative_swap_grid.back() - cumulative_swap_grid.front()) / ((double) cumulative_swap_grid.size() - 1);

            for (size_t ii = 1; ii < cumulative_swap_grid.size() - 1; ii++) {
                cumulative_swap_grid[ii] = cumulative_swap_grid[ii - 1] + step;
            }


            std::vector<double> new_temp_gradient(temp_gradient.size());
            new_temp_gradient.front() = temp_gradient.back();
            new_temp_gradient.back()  = 1.0;
            for (size_t ii = 1; ii < temp_gradient.size() - 1; ii++) {
                new_temp_gradient[ii] = spl.solve_one(cumulative_swap_grid[ii]);
            }

            std::reverse(new_temp_gradient.begin(), new_temp_gradient.end());

            for (size_t ii = 1; ii < temp_gradient.size() - 1; ii++) {
                chains[swap_indices[ii]].target->setTemperature(new_temp_gradient[ii]);
            }

            temp_gradient = new_temp_gradient;
            num_swaps     = 0;
            std::fill(swap_barriers.begin(), swap_barriers.end(), 0.0);
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

        Likelihood hotValueThreadSafe() {
            return chains[swap_indices[0]].target->valueThreadSafe();
        }

        void printModelLlik() {
            fmt::print("Current P/L: {0:.2f} / {1:.2f}. Post: {2:.2f}", chains[swap_indices[0]].model->getPrior(), chains[swap_indices[0]].model->getLikelihood(), chains[swap_indices[0]].target->value());
        }


        void finalize() {
            chains[swap_indices[0]].stateLogger->finalize();
            chains[swap_indices[0]].modelLogger->finalize();
        }


        std::vector<Chain> chains{};
        std::vector<int> swap_acceptance_rates{};
        std::vector<size_t> swap_indices{};

        std::vector<double> swap_barriers{};
        std::vector<double> temp_gradient{};
        std::vector<int> swap_store{};

        int num_swaps  = 0;
        bool even_swap = false;

        std::shared_ptr<boost::random::mt19937> r_;
        unsigned int num_cores_;
        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<unsigned int> seed_dist_{};
    };


}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_REPLICAEXCHANGE_H
