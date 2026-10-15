//
// Capturing a Checkpoint of an Ensemble and restoring an Ensemble from one: the State's Quantities as the Checkpoint's
// ChainState and back, with the ensemble recording and restoring the rest itself (see
// core/samplers/meta/EnsembleRecord.h). Model caches are not saved; restore rebuilds them from the State. The
// Checkpoint itself, and its file, are Checkpoint.h's: nothing here is needed to read or write one.
//

#ifndef TRANSMISSION_NETWORKS_APP_ENSEMBLECHECKPOINT_H
#define TRANSMISSION_NETWORKS_APP_ENSEMBLECHECKPOINT_H

#include "Checkpoint.h"
#include "Ensemble.h"
#include "State.h"

#include <nlohmann/json.hpp>

#include <memory>
#include <vector>

namespace transmission_nets::impl::Model::checkpoint {
    // The State's Quantities as a ChainState, and back into a State no Model has been built on yet, from the same input
    // as `layout`. restoreState throws Refused for a ChainState that doesn't fit the State.
    Layout stateLayout(const State& state);
    ChainState captureState(const State& state);
    void restoreState(const Layout& layout, const ChainState& saved, State& state);

    // Takes a checkpoint of `ensemble`.
    Checkpoint capture(const Ensemble& ensemble, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds);

    // Builds an ensemble from `checkpoint` and the Run's input, with `runRng` as its top-level engine, continued from
    // where the checkpoint left it, pending ladder barriers included; a new Run started from it calls
    // restartBurnIn(). Throws Refused for a checkpoint that doesn't fit the input. Whether it fits the Run is the Run
    // recipe's to check.
    std::unique_ptr<Ensemble> restore(const Checkpoint& checkpoint,
                                      const std::shared_ptr<Engine>& runRng,
                                      const nlohmann::json& input,
                                      const std::vector<core::computation::Probability>& symptomaticIDPDist,
                                      const std::vector<core::computation::Probability>& asymptomaticIDPDist,
                                      bool nullModel,
                                      unsigned int numCores);
}// namespace transmission_nets::impl::Model::checkpoint

#endif//TRANSMISSION_NETWORKS_APP_ENSEMBLECHECKPOINT_H
