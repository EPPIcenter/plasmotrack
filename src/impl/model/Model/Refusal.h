//
// A Run refused to start or resume: one kind per rule, and a message for the console. Every refusal exits 1.
//

#ifndef TRANSMISSION_NETWORKS_APP_REFUSAL_H
#define TRANSMISSION_NETWORKS_APP_REFUSAL_H

#include <stdexcept>
#include <string>
#include <utility>

namespace transmission_nets::impl::Model {

    enum class RefusalKind {
        // Which Run is meant.
        ResumeWithStartFrom,  // --resume and --start-from given together
        ResumeWithOutputDir,  // --resume given with an output directory
        NoOutputDir,          // a new Run given no output directory
        OutputHoldsRun,       // a new Run's output directory already holds a Run
        NothingToResume,      // the store given to --resume holds no Run
        // What was given, on the file system.
        InputMissing,
        InvalidInput,// the input does not describe a study
        OutputDirMissing,
        UnreadableRunManifest,
        UnusableCheckpoint,// can't be opened, isn't a checkpoint, has another format version, or doesn't fit the input
        // Resuming.
        RecipeChanged,  // an option given again differs from the Run's
        SamplesShrunk,  // -s below what the Run plans
        RunIdMismatch,  // the store's checkpoint is of another Run
        IdpChanged,     // an IDP differs from the Run's
        // A checkpoint resumed or started from.
        InputChanged,
        EngineVersionChanged,
        ChainCountChanged,
        // Value ranges.
        ThinOutOfRange,
        WindowOutOfRange,
        ParentSetFloorOutOfRange,
    };

    struct Refusal {
        RefusalKind kind;
        std::string message;
    };

    inline Refusal refuse(const RefusalKind kind, std::string message) {
        return {.kind = kind, .message = std::move(message)};
    }

    // Thrown where a refusal is found below the Run recipe: reading a stored Run or a checkpoint.
    class Refused : public std::runtime_error {
    public:
        explicit Refused(Refusal refusal) : std::runtime_error(refusal.message), refusal(std::move(refusal)) {}
        const Refusal refusal;
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_REFUSAL_H
