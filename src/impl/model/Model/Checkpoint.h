//
// The Checkpoint: everything a Run carries between Sampling iterations, so that it can be resumed or a new Run started
// from its State, and the codec that writes it as CBOR and reads it back. The ensemble's part is the ensemble's own
// record; the State's part is a ChainState of its Quantities. Taking one of an Ensemble and restoring an Ensemble from
// one are EnsembleCheckpoint.h's; this header needs no Model.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHECKPOINT_H
#define TRANSMISSION_NETWORKS_APP_CHECKPOINT_H

#include "impl/model/Model/ParentSetIds.h"
#include "impl/model/Model/Refusal.h"
#include "impl/model/Model/RunPosition.h"

#include "core/samplers/meta/EnsembleRecord.h"

#include <boost/random/mersenne_twister.hpp>

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
    // 5: the exchange state carries its burn-in iteration count.
    inline constexpr int FORMAT_VERSION = 5;

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

    // The ensemble's own record, with each chain's State as a ChainState; see core/samplers/meta/EnsembleRecord.h.
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

    // A checkpoint that can't be used is a refusal, which exits 1.
    [[noreturn]] void refuse(std::string message);

    // Writes `checkpoint` to `file` as CBOR, replacing any file there durably (see core/io/DurableFile.h). A failure
    // ends the process with OUTPUT_FAILURE, naming the path.
    void save(const Checkpoint& checkpoint, const std::filesystem::path& file);

    // Throws Refused for a file that can't be opened, isn't a checkpoint or has another format version.
    Checkpoint load(const std::filesystem::path& file);

}// namespace transmission_nets::impl::Model::checkpoint

#endif//TRANSMISSION_NETWORKS_APP_CHECKPOINT_H
