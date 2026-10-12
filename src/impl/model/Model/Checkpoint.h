//
// Saves and restores everything a Run carries between Sampling iterations, so that a Run can be resumed or a new Run
// started from a saved State. Written as CBOR. Model caches are not saved; restore rebuilds them from the State.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINT_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINT_H

#include "core/samplers/meta/ReplicaExchange.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/ParentSetIds.h"
#include "impl/model/Model/Refusal.h"
#include "impl/model/Model/RunPosition.h"
// SequentialScheduler.h relies on the sampler headers SampleScheduler.h pulls in.
#include "impl/model/Model/SampleScheduler.h"
#include "impl/model/Model/SequentialScheduler.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {
    using Ensemble = core::samplers::ReplicaExchange<State, Model, SequentialSampleScheduler>;
}

namespace transmission_nets::impl::Model::checkpoint {
    using Engine = boost::random::mt19937;

    // Written into every checkpoint. Raise it when a change means an older engine would misread the file.
    inline constexpr int FORMAT_VERSION = 2;

    // The checkpoints in a Run's store: the latest Window's, and the one at burn-in end.
    inline constexpr const char* LATEST_FILE  = "checkpoint.cbor";
    inline constexpr const char* BURN_IN_FILE = "checkpoint-burn-in.cbor";

    // Which Run a checkpoint belongs to, and what it was made from.
    struct Identity {
        std::string runId;
        std::string inputHash;
        std::string engineVersion;

        bool operator==(const Identity&) const = default;
    };

    // Where the Run stood when the checkpoint was taken. A checkpoint is taken only once a Window is committed, so it
    // always has a last closed Window.
    using transmission_nets::impl::Model::Phase;
    using transmission_nets::impl::Model::RunPosition;

    // The Parent-set id dictionary the Run recorder owns and the checkpoint carries; see ParentSetIds.h.
    using transmission_nets::impl::Model::ParentSetIds;
    using transmission_nets::impl::Model::ParentSetIdDictionary;

    // What every chain shares, named once per checkpoint.
    struct Layout {
        std::vector<std::string> infections;   // Infection ids, in input order
        std::vector<std::string> latentParents;// Latent parent ids, in input order
        std::vector<std::string> loci;         // Locus labels, sorted
        std::vector<std::string> samplers;     // sampler ids, in the order the schedulers registered them

        bool operator==(const Layout&) const = default;
    };

    // One chain's full State, in Layout order. Genotypes are packed Allele bits: for each Infection (or Latent parent),
    // for each of its Loci, one bit per Allele, filling each byte from its most significant bit.
    struct ChainState {
        std::vector<std::uint8_t> latentGenotypes;
        std::vector<std::uint8_t> latentParentGenotypes;
        std::vector<std::vector<double>> alleleFrequencies;// by Locus
        std::vector<std::uint32_t> ordering;               // positions in Layout::infections, first acquired first
        std::vector<double> infectionDurations;            // by Infection
        std::vector<double> expectedFalsePositives;        // by Infection
        std::vector<double> expectedFalseNegatives;        // by Infection
        double meanCoi = 0;
        double meanStrainsTransmitted = 0;
        double parentSetSizeProb = 0;

        bool operator==(const ChainState&) const = default;
    };

    // One replica: its State, its RNG engine's state, its samplers' counts and adaptation state, and how many Sampling
    // iterations its scheduler has run, which decides when each sampler adapts.
    struct Replica {
        ChainState state;
        std::string rng;
        std::vector<core::samplers::SamplerCounts> samplerCounts;              // by sampler, in Layout order
        std::map<std::uint32_t, core::samplers::AdaptationState> adaptation;   // by position in Layout::samplers; only samplers that adapt
        std::int64_t schedulerIterations = 0;

        bool operator==(const Replica&) const = default;
    };

    struct Checkpoint {
        Identity identity;
        RunPosition position;
        ParentSetIdDictionary parentSetIds;
        Layout layout;
        std::string runRng;               // the top-level engine, which draws replica exchange swaps
        std::vector<Replica> replicas;    // by replica, not by temperature slot
        Ensemble::ExchangeState exchange; // the temperature slots: replica, temperature, swap barriers and counts; its
                                          // burn-in iteration count isn't saved, and restore takes it from `position`

        bool operator==(const Checkpoint&) const = default;
    };

    // Takes a checkpoint of `ensemble`, whose top-level RNG engine is `runRng`.
    Checkpoint capture(const Ensemble& ensemble, const Engine& runRng, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds);

    // Writes `checkpoint` to `file` as CBOR, replacing any file there durably (see core/io/DurableFile.h). A failure
    // ends the process with OUTPUT_FAILURE, naming the path.
    void save(const Checkpoint& checkpoint, const std::filesystem::path& file);

    // Throws Refused for a file that can't be opened, isn't a checkpoint or has another format version.
    Checkpoint load(const std::filesystem::path& file);

    // Builds an ensemble from `checkpoint` and the Run's input, and sets `runRng` to the saved top-level engine. Throws
    // Refused for a checkpoint that doesn't fit the input. Whether it fits the Run is the Run recipe's to check.
    std::unique_ptr<Ensemble> restore(const Checkpoint& checkpoint,
                                      const std::shared_ptr<Engine>& runRng,
                                      const nlohmann::json& input,
                                      const std::vector<core::computation::Probability>& symptomaticIDPDist,
                                      const std::vector<core::computation::Probability>& asymptomaticIDPDist,
                                      bool nullModel,
                                      unsigned int numCores);
}// namespace transmission_nets::impl::Model::checkpoint

#endif//TRANSMISSION_NETWORKS_APP_CHECKPOINT_H
