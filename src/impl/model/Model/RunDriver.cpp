//
// The Run driver. See RunDriver.h.
//

#include "RunDriver.h"

namespace transmission_nets::impl::Model {
    RunDriver::RunDriver(Ensemble& ensemble, RunRecorder& recorder, std::filesystem::path store, checkpoint::Identity identity, Progress& progress, RunRecorder::Clock clock)
        : ensemble_(ensemble), recorder_(recorder), store_(std::move(store)), identity_(std::move(identity)), progress_(progress), clock_(std::move(clock)) {}

    RunOutcome RunDriver::drive(const RunStart& start, const RunPlan& plan, const StopQuery& stop) {
        // Every call to the recorder holds the lock, so a forced end never comes while one is in progress.
        const auto recording = [&](auto&& call) {
            const std::lock_guard held(lock_);
            return call();
        };

        // Sampling iterations already run in each Phase: none unless resuming.
        const auto* resumed             = std::get_if<ResumedRun>(&start);
        const bool resumedInSampling    = resumed and resumed->position.phase == Phase::Sampling;
        const std::int64_t burnInDone   = resumed ? (resumedInSampling ? plan.burnIn : resumed->position.iteration) : 0;
        const std::int64_t samplingDone = resumedInSampling ? resumed->position.iteration - plan.burnIn : 0;

        // Where the Run stands as of the last Window committed, or where it started.
        RunPosition at = resumed ? resumed->position : RunPosition{.lastClosedWindow = -1, .burnInRemaining = plan.burnIn, .samplingRemaining = plan.samples};
        if (resumed) {
            // The plan may have grown since the Checkpoint was taken.
            at.samplingRemaining = plan.samples - samplingDone;
        }

        // When progress was last reported with a line; a heartbeat goes out once HEARTBEAT_INTERVAL has passed since.
        const auto startedAt = clock_();
        auto lastReport      = startedAt;
        std::int64_t iterationsRun = 0;// Sampling iterations run by this call, in either Phase

        // The Checkpoint of the complete in-flight Window, taken before the Sampling iteration that commits it.
        std::optional<checkpoint::Checkpoint> windowEnd;
        recording([&] {
            driving_       = true;
            firstRecorded_ = false;
            recorder_.onWindowCommitted([&](const WindowCommit& commit) {
                auto taken = windowEnd ? std::move(*windowEnd) : checkpoint::capture(ensemble_, identity_, {}, {});
                windowEnd.reset();
                taken.position     = commit.position;
                taken.parentSetIds = commit.parentSetIds;
                checkpoint::save(taken, store_ / checkpoint::LATEST_FILE);
                if (commit.burnInEnded) {
                    checkpoint::save(taken, store_ / checkpoint::BURN_IN_FILE);
                }
                at = commit.position;
                progress_.windowClosed(at, commit.summary);
                lastReport = clock_();
            });
        });
        // However drive() is left, a later forced end records nothing and the recorder calls back into nothing here.
        struct Finish {
            RunDriver& driver;
            ~Finish() {
                const std::lock_guard held(driver.lock_);
                driver.driving_ = false;
                driver.recorder_.onWindowCommitted({});
            }
        } const finish{*this};

        progress_.started(at, resumed != nullptr);
        recording([&] {
            if (resumed) {
                recorder_.recordResume({.position     = resumed->position,
                                        .checkpoint   = (store_ / checkpoint::LATEST_FILE).string(),
                                        .samples      = plan.samples,
                                        .parentSetIds = resumed->parentSetIds},
                                       ensemble_.snapshot());
            } else {
                recorder_.recordStart(ensemble_.snapshot());
            }
        });

        // Where the Run stands once it has run Sampling iteration `iteration` of the Run, in `phase`.
        const auto after = [&](const Phase phase, const std::int64_t iteration) {
            const bool burningIn = phase == Phase::BurnIn;
            return RunPosition{.iteration         = iteration,
                               .phase             = phase,
                               .lastClosedWindow  = at.lastClosedWindow,
                               .burnInRemaining   = burningIn ? plan.burnIn - iteration : 0,
                               .samplingRemaining = burningIn ? plan.samples : plan.samples - (iteration - plan.burnIn)};
        };

        // The one Sampling iteration body, for both Phases. `iteration` counts within the Run, from 1.
        const auto iterate = [&](const Phase phase, const std::int64_t iteration) {
            recording([&] {
                if (recorder_.windowComplete()) {
                    windowEnd = checkpoint::capture(ensemble_, identity_, {}, {});
                }
            });
            std::optional<TemperatureAdaptation> adaptation;
            if (phase == Phase::BurnIn) {
                adaptation = ensemble_.burnin();
            } else {
                ensemble_.sample();
            }
            const auto snapshot = ensemble_.snapshot();
            recording([&] { recorder_.recordIteration(snapshot); });
            if (adaptation) {
                recording([&] { recorder_.recordTemperatureAdaptation(*adaptation); });
                // Not a report that puts the heartbeat off: the console prints no line for it.
                progress_.ladderAdapted(after(phase, iteration), *adaptation);
            }
            ++iterationsRun;
            if (const auto now = clock_(); now - lastReport >= HEARTBEAT_INTERVAL) {
                const auto& hot = snapshot.slots.front();
                progress_.heartbeat(after(phase, iteration),
                                    {.hotPrior = hot.prior, .hotLikelihood = hot.likelihood, .samplingIterations = iterationsRun, .seconds = std::chrono::duration<double>(now - startedAt).count()});
                lastReport = now;
            }
        };

        // A stop request takes effect before the next Sampling iteration.
        std::string signal;
        const auto stopRequested = [&] {
            const auto asked = stop();
            if (!asked) {
                return false;
            }
            signal = *asked;
            recording([&] {
                recorder_.recordInterrupt(signal, Interrupt::First);
                firstRecorded_ = true;
            });
            return true;
        };

        bool stopped = false;
        for (std::int64_t kk = burnInDone; kk < plan.burnIn and !(stopped = stopRequested()); ++kk) {
            iterate(Phase::BurnIn, kk + 1);
        }
        if (!stopped and !resumedInSampling) {
            // Adaptation stops here. The burn-in checkpoint is taken when the Window holding burn-in end commits.
            recording([&] { recorder_.recordBurnInEnd(ensemble_.finalAdaptation()); });
        }
        for (std::int64_t jj = samplingDone; jj < plan.samples and !stopped and !(stopped = stopRequested()); ++jj) {
            iterate(Phase::Sampling, plan.burnIn + jj + 1);
        }
        recording([&] { recorder_.recordEnd(); });

        if (stopped) {
            progress_.stopped(at, signal);
            return Stopped{signal};
        }
        progress_.ended(at);
        return Completed{};
    }

    void RunDriver::forceEnd(const std::string& firstSignal, const std::string& signal) {
        // Never released: the process ends before anything else reaches the recorder.
        lock_.lock();
        if (!driving_) {
            return;
        }
        if (!firstRecorded_) {
            recorder_.recordInterrupt(firstSignal, Interrupt::First);
        }
        recorder_.recordInterrupt(signal, Interrupt::Second);
    }
}// namespace transmission_nets::impl::Model
