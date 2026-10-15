//
// Capturing a Checkpoint of an Ensemble and restoring an Ensemble from one. See EnsembleCheckpoint.h.
//

#include "EnsembleCheckpoint.h"

#include <fmt/core.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace transmission_nets::impl::Model::checkpoint {
    namespace {
        // The saved value for `key`; a missing key means the checkpoint was made from another input.
        template<typename Value>
        const Value& savedFor(const std::map<std::string, Value>& saved, const std::string& key, const std::string& what) {
            const auto it = saved.find(key);
            if (it == saved.end()) {
                refuse(fmt::format("checkpoint has no {} {}; it was made from another input", what, key));
            }
            return it->second;
        }

        void requireSameNames(const std::vector<std::string>& saved, const std::vector<std::string>& current, const std::string& what) {
            if (saved != current) {
                refuse(fmt::format("checkpoint has other {} than this input; it was made from another input", what));
            }
        }

        template<typename Value>
        void requireLength(const std::vector<Value>& saved, const size_t expected, const std::string& what) {
            if (saved.size() != expected) {
                refuse(fmt::format("checkpoint has {} {}, not {}; it was made from another input", saved.size(), what, expected));
            }
        }

        std::vector<std::string> infectionIds(const std::vector<std::shared_ptr<InfectionEvent>>& infections) {
            std::vector<std::string> ids;
            ids.reserve(infections.size());
            for (const auto& infection : infections) {
                ids.push_back(infection->id());
            }
            return ids;
        }

        std::vector<std::string> locusLabels(const State& state) {
            std::vector<std::string> labels;
            for (const auto& [label, locus] : state.loci) {
                labels.push_back(label);
            }
            return labels;
        }

        // Appends Allele bits to bytes, filling each byte from its most significant bit.
        class BitWriter {
        public:
            void push(const bool bit) {
                if (bits_ % 8 == 0) {
                    bytes_.push_back(0);
                }
                bytes_.back() |= static_cast<std::uint8_t>(static_cast<unsigned>(bit) << (7 - bits_ % 8));
                ++bits_;
            }
            std::vector<std::uint8_t> take() { return std::move(bytes_); }

        private:
            std::vector<std::uint8_t> bytes_;
            size_t bits_ = 0;
        };

        class BitReader {
        public:
            explicit BitReader(const std::vector<std::uint8_t>& bytes) : bytes_(bytes) {}
            bool next() {
                if (bits_ / 8 >= bytes_.size()) {
                    refuse("checkpoint has fewer Allele bits than this input's genotypes; it was made from another input");
                }
                const bool bit = (bytes_[bits_ / 8] >> (7 - bits_ % 8)) & 1;
                ++bits_;
                return bit;
            }
            void requireAllRead() const {
                if ((bits_ + 7) / 8 != bytes_.size()) {
                    refuse("checkpoint has more Allele bits than this input's genotypes; it was made from another input");
                }
            }

        private:
            const std::vector<std::uint8_t>& bytes_;
            size_t bits_ = 0;
        };

        // Visits a State's Quantities into a ChainState.
        struct Capturing {
            const State& state;
            ChainState saved{};
            BitWriter genotypes{};

            void scalar(const char* name, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                saved.scalars.emplace(name, p->value());
            }
            void perInfection(const char* name, size_t, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                saved.perInfection[name].push_back(p->value());
            }
            void alleleFrequencies(const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<AlleleFrequencyImpl>>& p) {
                saved.alleleFrequencies.push_back(p->value().frequencies());
            }
            void latentGenotype(size_t, const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<Genotype>>& p) {
                const auto& alleles = p->value();
                for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                    genotypes.push(alleles.allele(allele));
                }
            }
            void ordering(const std::shared_ptr<OrderingImpl>& ordering) {
                std::map<const InfectionEvent*, std::uint32_t> position;
                for (std::uint32_t ii = 0; ii < state.infections.size(); ++ii) {
                    position.emplace(state.infections[ii].get(), ii);
                }
                for (const auto& infection : ordering->value()) {
                    saved.ordering.push_back(position.at(infection.get()));
                }
            }
            ChainState take() {
                saved.latentGenotypes = genotypes.take();
                return std::move(saved);
            }
        };

        // Visits a ChainState into a State no Model has been built on yet.
        struct Restoring {
            State& state;
            const ChainState& saved;
            BitReader genotypes;
            size_t locusAt = 0;

            void scalar(const char* name, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                p->initializeValue(savedFor(saved.scalars, name, "Quantity"));
            }
            void perInfection(const char* name, const size_t index, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                const auto& values = savedFor(saved.perInfection, name, "per-Infection Quantity");
                if (index == 0) {
                    requireLength(values, state.infections.size(), name);
                }
                p->initializeValue(values[index]);
            }
            void alleleFrequencies(const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<AlleleFrequencyImpl>>& p) {
                p->initializeValue(AlleleFrequencyImpl::exactly(saved.alleleFrequencies[locusAt++]));
            }
            void latentGenotype(size_t, const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<Genotype>>& p) {
                Genotype alleles(static_cast<int>(p->value().totalAlleles()));
                for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                    alleles.set(allele, genotypes.next());
                }
                p->initializeValue(alleles);
            }
            void ordering(const std::shared_ptr<OrderingImpl>& ordering) {
                genotypes.requireAllRead();
                std::vector<std::shared_ptr<InfectionEvent>> order;
                for (const auto at : saved.ordering) {
                    if (at >= state.infections.size()) {
                        refuse(fmt::format("checkpoint's Ordering has Infection {} of {}", at, state.infections.size()));
                    }
                    order.push_back(state.infections[at]);
                }
                ordering->restoreOrder(std::move(order));
            }
        };

    }// namespace

    Layout stateLayout(const State& state) {
        return {.infections = infectionIds(state.infections), .latentParents = infectionIds(state.latentParents), .loci = locusLabels(state)};
    }

    ChainState captureState(const State& state) {
        Capturing capturing{.state = state};
        state.quantities(capturing);
        return capturing.take();
    }

    // Fills a State no Model has been built on yet, from the same input as `layout`.
    void restoreState(const Layout& layout, const ChainState& saved, State& state) {
        requireSameNames(layout.infections, infectionIds(state.infections), "Infections");
        requireSameNames(layout.latentParents, infectionIds(state.latentParents), "Latent parents");
        requireSameNames(layout.loci, locusLabels(state), "Loci");
        requireLength(saved.alleleFrequencies, state.loci.size(), "Allele frequencies");
        requireLength(saved.ordering, state.infections.size(), "Infections in the Ordering");
        Restoring restoring{.state = state, .saved = saved, .genotypes = BitReader(saved.latentGenotypes)};
        state.quantities(restoring);
    }


    Checkpoint capture(const Ensemble& ensemble, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds) {
        Checkpoint checkpoint{.identity = std::move(identity), .position = position, .parentSetIds = std::move(parentSetIds), .layout = stateLayout(*ensemble.chains.front().state)};
        checkpoint.ensemble = ensemble.record([&](const State& state) {
            if (infectionIds(state.infections) != checkpoint.layout.infections) {
                throw std::logic_error("chains of one ensemble have different Infections");
            }
            return captureState(state);
        });
        return checkpoint;
    }

    std::unique_ptr<Ensemble> restore(const Checkpoint& checkpoint, const std::shared_ptr<Engine>& runRng, const Study& study, const bool nullModel, const unsigned int numCores) {
        std::unique_ptr<Ensemble> ensemble;
        try {
            ensemble = std::make_unique<Ensemble>(checkpoint.ensemble, runRng, numCores, [&](const size_t, const ChainState& saved, std::shared_ptr<Engine> chainRng) {
                auto state = std::make_shared<State>(study, std::move(chainRng), nullModel);
                restoreState(checkpoint.layout, saved, *state);
                return state;
            });
        } catch (const std::invalid_argument& e) {
            refuse(fmt::format("checkpoint of Run {} does not fit this input: {}", checkpoint.identity.runId, e.what()));
        }

        return ensemble;
    }
}// namespace transmission_nets::impl::Model::checkpoint
