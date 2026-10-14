//
// A Proposal: one Metropolis-Hastings step of a sampler, from the first parameter it changes to its accept or reject.
// The parameters it touches, the target's value before and after, and whether it was accepted live here; a sampler
// decides what to change and the Hastings adjustment, nothing else.
//
// The rules every step must follow live here rather than in the samplers:
// - a step's saved states carry an id no other open step shares, so steps nest (a multiple-try move evaluates tries
//   while its own change is saved) and a restore pops exactly its own save;
// - the target is evaluated before a change is accepted or restored: an accept marks the model's terms clean, and a
//   term may record a change when it happens and apply it when it is next read;
// - a change inside a nested step leaves the enclosing steps unevaluated, so they evaluate again before they settle;
// - the accept/reject decision draws its uniform after the proposed state is evaluated;
// - every step counts once in the sampler's counts, accepted or not, unless it is a try that is reverted.
//
// A sampler's target must offer value(); isDirty() is checked in debug builds when it exists.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROPOSAL_H
#define TRANSMISSION_NETWORKS_APP_PROPOSAL_H

#include "core/samplers/SamplerCounts.h"

#include <boost/container/small_vector.hpp>
#include <boost/random/uniform_01.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <utility>

namespace transmission_nets::core::samplers {

    namespace detail {
        // The Proposals open on this thread, innermost last. Ids come from here too.
        struct OpenProposals {
            boost::container::small_vector<bool*, 4> evaluated{};
            std::uint32_t nextId = 1;

            static OpenProposals& here() {
                static thread_local OpenProposals open;
                return open;
            }
            // Positive and distinct from every id another open Proposal on this thread can hold.
            int takeId() {
                return static_cast<int>(nextId++ % 0x7ffffffeu) + 1;
            }
        };
    }// namespace detail

    template<typename Target, typename Engine>
    class Proposal {
    public:
        // Reads the target at the current state. `counts` is the sampler's; accept() and reject() count into it.
        Proposal(Target& target, Engine& rng, SamplerCounts& counts) : target_(target), rng_(rng), counts_(counts), id_(open().takeId()), current_(target.value()) {
            open().evaluated.push_back(&evaluated_);
        }
        Proposal(const Proposal&)            = delete;
        Proposal& operator=(const Proposal&) = delete;
        ~Proposal() {
            assert(settled_ && "a Proposal ends in accept(), reject() or revert()");
            open().evaluated.pop_back();
        }

        // The target's value before any change.
        [[nodiscard]] double current() const noexcept {
            return current_;
        }

        // Saves `parameter` under this Proposal's id and sets it. A parameter may be changed more than once.
        template<typename Parameter, typename Value>
        void change(Parameter& parameter, Value&& value) {
            parameter.saveState(id_);
            touched_.push_back({&parameter,
                                [](void* p) { static_cast<Parameter*>(p)->acceptState(); },
                                [](void* p, const int id) { static_cast<Parameter*>(p)->restoreState(id); }});
            parameter.setValue(std::forward<Value>(value));
            for (bool* evaluated : open().evaluated) {
                *evaluated = false;
            }
        }

        // The target's value at the proposed state.
        double value() {
            const double proposed = target_.value();
            if constexpr (requires { target_.isDirty(); }) {
                assert(!target_.isDirty());
            }
            evaluated_ = true;
            return proposed;
        }

        // Keeps the change and counts an accepted proposal. Not inside another open Proposal: an accept clears every
        // save a parameter holds, the outer step's included.
        void accept() {
            assert(open().evaluated.size() == 1 && "accept() inside an open outer Proposal would commit its change too");
            evaluateIfNeeded();
            for (const auto& touched : touched_) {
                touched.accept(touched.parameter);
            }
            settle(true);
        }

        // Restores every parameter changed and counts a rejected proposal.
        void reject() {
            restoreAll();
            settle(false);
        }

        // Restores every parameter changed without counting: a try evaluated on the way to a move.
        void revert() {
            restoreAll();
            settled_ = true;
        }

        // Draws the accept/reject uniform, then accepts or rejects. Returns whether it accepted.
        bool decide(const double logAcceptProbability) {
            boost::random::uniform_01<> uniform;
            const bool accepted = std::log(uniform(rng_)) <= logAcceptProbability;
            if (accepted) {
                accept();
            } else {
                reject();
            }
            return accepted;
        }

        // Metropolis-Hastings: evaluates the proposed state and decides on its log ratio to the current one plus the
        // adjustment.
        bool metropolis(const double logHastingsAdjustment = 0) {
            const double proposed = value();
            return decide(proposed - current_ + logHastingsAdjustment);
        }

    private:
        struct Touched {
            void* parameter;
            void (*accept)(void*);
            void (*restore)(void*, int);
        };

        static detail::OpenProposals& open() {
            return detail::OpenProposals::here();
        }

        void evaluateIfNeeded() {
            if (!evaluated_) {
                value();
            }
        }

        void restoreAll() {
            evaluateIfNeeded();
            for (auto it = touched_.rbegin(); it != touched_.rend(); ++it) {
                it->restore(it->parameter, id_);
            }
        }

        void settle(const bool accepted) noexcept {
            ++counts_.proposed;
            counts_.accepted += accepted;
            settled_ = true;
        }

        Target& target_;
        Engine& rng_;
        SamplerCounts& counts_;
        int id_;
        double current_;
        bool evaluated_ = true;
        bool settled_   = false;
        boost::container::small_vector<Touched, 8> touched_{};
    };

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_PROPOSAL_H
