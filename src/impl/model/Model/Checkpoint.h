//
// Saves and restores everything a Run carries between Sampling iterations, so that a Run can be resumed or a new Run
// started from a saved State. Written as CBOR. Model caches are not saved; restore rebuilds them from the State.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINT_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINT_H

#include "impl/model/Model/Ensemble.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/ParentSetIds.h"
#include "impl/model/Model/Refusal.h"
#include "impl/model/Model/RunPosition.h"
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

namespace transmission_nets::impl::Model::checkpoint {
    using Engine = boost::random::mt19937;

    // Written into every checkpoint. Raise it when a change means an older engine would misread the file.
    // 3: latent parents are summed out of the model, so their genotypes are no longer saved.
    // 4: Quantities are saved by name: the scalars and the per-Infection values are maps keyed by Quantity name.
    inline constexpr int FORMAT_VERSION = 4;

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

    // What every chain's State shares, named once per checkpoint. The sampler ids every chain shares are the ensemble
    // record's.
    struct Layout {
        std::vector<std::string> infections;   // Infection ids, in input order
        std::vector<std::string> latentParents;// Latent parent ids, in input order
        std::vector<std::string> loci;         // Locus labels, sorted

        bool operator==(const Layout&) const = default;
    };

    // One chain's Quantities (see State::quantities), in Layout order. Genotypes are packed Allele bits: for each
    // Infection, for each of its Loci, one bit per Allele, filling each byte from its most significant bit.
    struct ChainState {
        std::map<std::string, double> scalars;                  // by Quantity name
        std::map<std::string, std::vector<double>> perInfection;// by Quantity name, then by Infection
        std::vector<std::vector<double>> alleleFrequencies;     // by Locus
        std::vector<std::uint8_t> latentGenotypes;
        std::vector<std::uint32_t> ordering;// positions in Layout::infections, first acquired first
        bool operator==(const ChainState&) const = default;
    };

    // The ensemble's own record, with each chain's State as a ChainState; see core/samplers/meta/EnsembleRecord.h. Its
    // exchange state's burn-in iteration count isn't saved: restore takes it from the Run position.
    using EnsembleRecord = core::samplers::EnsembleRecord<ChainState>;
    using Replica        = EnsembleRecord::Replica;

    struct Checkpoint {
        Identity identity;
        RunPosition position;
        ParentSetIdDictionary parentSetIds;
        Layout layout;
        EnsembleRecord ensemble;

        bool operator==(const Checkpoint&) const = default;
    };

    // The State's Quantities as a ChainState, and back into a State no Model has been built on yet, from the same input
    // as `layout`. restoreState throws Refused for a ChainState that doesn't fit the State.
    Layout stateLayout(const State& state);
    ChainState captureState(const State& state);
    void restoreState(const Layout& layout, const ChainState& saved, State& state);

    // Takes a checkpoint of `ensemble`.
    Checkpoint capture(const Ensemble& ensemble, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds);

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
