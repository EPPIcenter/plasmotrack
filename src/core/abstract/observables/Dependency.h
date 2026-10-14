//
// Dependency wiring: the one place that says what it means for a term or a cached computation to depend on an input.
//
// A dependent of an input is marked dirty when the input changes, after any bookkeeping of its own, and follows the
// input's saved states: it saves with it, accepts with it and restores with it. Its restored value is the one it had
// when the input was saved, so a restore leaves it clean; a save evaluates it first (see Checkpointable::saveState), so
// that value is always current. An accept settles whatever the dependent was evaluated to and marks it clean; the
// Proposal evaluates the target before every accept, so a dependent reached through the Model is never accepted unread.
//
// An input is a Parameter, a container that passes its Parameters' events through, or a cached computation: the first
// two announce post_change, the last set_dirty. The compiler picks.
//

#ifndef TRANSMISSION_NETWORKS_APP_DEPENDENCY_H
#define TRANSMISSION_NETWORKS_APP_DEPENDENCY_H

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

namespace transmission_nets::core::abstract {
    namespace detail {
        template<typename T>
        struct IsSharedPtr : std::false_type {};
        template<typename T>
        struct IsSharedPtr<std::shared_ptr<T>> : std::true_type {};

        // The input itself, whether given directly or through a shared_ptr.
        template<typename Input>
        decltype(auto) deref(Input& input) {
            if constexpr (IsSharedPtr<std::remove_cvref_t<Input>>::value) {
                return *input;
            } else {
                return (input);
            }
        }

        template<typename Input>
        concept AnnouncesChanges = requires(Input& input, std::function<void()> cb) { input.add_post_change_listener(cb); };
        template<typename Input>
        concept AnnouncesDirt = requires(Input& input, std::function<void()> cb) { input.add_set_dirty_listener(cb); };
        template<typename Input>
        concept HasSavedStates = requires(Input& input, std::function<void(int)> onId, std::function<void()> on) {
            input.add_save_state_listener(onId);
            input.add_accept_state_listener(on);
            input.add_restore_state_listener(onId);
        };

        // `dependent` saves, accepts and restores with `input`, and is clean after an accept or a restore.
        template<typename Dependent, typename Input>
        void followSavedStates(Dependent* dependent, Input& input) {
            static_assert(HasSavedStates<Input>, "an input has saved states to follow");
            input.add_save_state_listener([dependent](const int savedStateId) { dependent->saveState(savedStateId); });
            input.add_accept_state_listener([dependent]() {
                dependent->acceptState();
                dependent->setClean();
            });
            input.add_restore_state_listener([dependent](const int savedStateId) {
                dependent->restoreState(savedStateId);
                dependent->setClean();
            });
        }
    }// namespace detail

    // `dependent` depends on `input`: when the input changes, `onChange` runs and the dependent is marked dirty, and it
    // follows the input's saved states.
    template<typename Dependent, typename InputLike, typename OnChange>
        requires std::invocable<OnChange&>
    void dependsOn(Dependent* dependent, InputLike&& inputLike, OnChange&& onChange) {
        auto& input      = detail::deref(inputLike);
        using Input      = std::remove_cvref_t<decltype(input)>;
        auto changed     = [dependent, onChange = std::forward<OnChange>(onChange)]() mutable {
            onChange();
            dependent->setDirty();
        };
        static_assert(detail::AnnouncesChanges<Input> != detail::AnnouncesDirt<Input>, "an input announces either changes or dirt");
        if constexpr (detail::AnnouncesDirt<Input>) {
            input.add_set_dirty_listener(changed);
        } else {
            input.add_post_change_listener(changed);
        }
        detail::followSavedStates(dependent, input);
    }

    // `dependent` depends on each of `inputs`, with nothing of its own to do when one changes.
    template<typename Dependent, typename... InputLikes>
        requires(sizeof...(InputLikes) > 0 and (not std::invocable<InputLikes&> and ...))
    void dependsOn(Dependent* dependent, InputLikes&&... inputs) {
        (dependsOn(dependent, inputs, []() {}), ...);
    }

    // `dependent` follows the saved states of `input` without being marked dirty by it: for a dependent told of the
    // input's changes some other way, or one whose cache records the changes itself.
    template<typename Dependent, typename InputLike>
    void savesWith(Dependent* dependent, InputLike&& inputLike) {
        detail::followSavedStates(dependent, detail::deref(inputLike));
    }

    // `container` passes `part`'s changes and saved states through as its own, so that a dependent of the container
    // depends on every part.
    template<typename Container, typename PartLike>
    void passThrough(Container* container, PartLike&& partLike) {
        auto& part = detail::deref(partLike);
        part.add_pre_change_listener([container]() { container->notify_pre_change(); });
        part.add_post_change_listener([container]() { container->notify_post_change(); });
        part.add_save_state_listener([container](const int savedStateId) { container->notify_save_state(savedStateId); });
        part.add_accept_state_listener([container]() { container->notify_accept_state(); });
        part.add_restore_state_listener([container](const int savedStateId) { container->notify_restore_state(savedStateId); });
    }
}// namespace transmission_nets::core::abstract

#endif//TRANSMISSION_NETWORKS_APP_DEPENDENCY_H
