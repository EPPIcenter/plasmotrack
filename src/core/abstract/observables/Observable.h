//
// Created by Maxwell Murphy on 1/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVABLE_H
#define TRANSMISSION_NETWORKS_APP_OBSERVABLE_H

#include <iostream>
#include <boost/container/flat_map.hpp>

#include "core/abstract/crtp.h"

//// Enables add and removing of callbacks for a given callback_type
#define CREATE_EVENT(callback_name, callback_type)                                          \
public:                                                                                     \
    template<typename ...Args>                                                              \
    void notify_##callback_name(Args... args) const noexcept {                              \
        this->notify(this->callback_name##_callbacks_, args...);                            \
    }                                                                                       \
                                                                                            \
    auto add_##callback_name##_listener(const callback_type& cb) noexcept -> ListenerId {   \
        return this->add_listener(cb, this->callback_name##_callbacks_);                    \
    }                                                                                       \
                                                                                            \
    auto remove_##callback_name##_listener(const ListenerId_t id) noexcept -> bool {          \
        return this->remove_listener(id, this->callback_name##_callbacks_);                 \
    }                                                                                       \
private:                                                                                    \
    ObserverMap<ListenerId_t, callback_type> callback_name##_callbacks_{};                  \

#define CREATE_KEYED_EVENT(callback_name, KeyType, CallbackType)                                                    \
public:                                                                                                             \
    template<typename ...Args>                                                                                      \
    void keyed_notify_##callback_name(KeyType key, Args... args) const noexcept {                                   \
        this->keyed_notify(key, this->keyed_##callback_name##_callbacks_, args...);                                 \
    }                                                                                                               \
                                                                                                                    \
    auto add_keyed_##callback_name##_listener(KeyType key, const CallbackType& cb) noexcept -> ListenerId {         \
        return this->add_keyed_listener(key, cb, this->keyed_##callback_name##_callbacks_);                         \
    }                                                                                                               \
                                                                                                                    \
    auto remove_keyed_##callback_name##_listener(KeyType key, const ListenerId_t id) noexcept -> bool {               \
        return this->remove_keyed_listener(key, id, this->keyed_##callback_name##_callbacks_);                      \
    }                                                                                                               \
                                                                                                                    \
    void register_##callback_name##_listener_key(KeyType key) noexcept {                                            \
        this->keyed_##callback_name##_callbacks_.emplace(key, ObserverMap<ListenerId_t, CallbackType>{});           \
    }                                                                                                               \
                                                                                                                    \
private:                                                                                                            \
    ObserverMap<KeyType, ObserverMap<ListenerId_t, CallbackType>> keyed_##callback_name##_callbacks_{};             \

#define CRTP_CREATE_EVENT(callback_name, callback_type)                                     \
public:                                                                                     \
    void notify_##callback_name() const noexcept {                                          \
        this->underlying().notify(this->callback_name##_callbacks_);                        \
    }                                                                                       \
                                                                                            \
    auto add_##callback_name##_listener(const callback_type& cb) noexcept -> ListenerId {   \
        return this->underlying().add_listener(cb, this->callback_name##_callbacks_);       \
    }                                                                                       \
                                                                                            \
    auto remove_##callback_name##_listener(const ListenerId_t id) noexcept -> bool {        \
        return this->underlying().remove_listener(id, this->callback_name##_callbacks_);    \
    }                                                                                       \
protected:                                                                                  \
    ObserverMap<ListenerId_t, callback_type> callback_name##_callbacks_{};                  \


// Data structures to support listener storage
template<typename Key, typename Value> using ObserverMap = boost::container::flat_map<Key, Value>;

using ListenerId_t = uint_fast32_t;
enum ListenerId : ListenerId_t {
};

template<typename T>
class Observable : public crtp<T, Observable> {
public:
    static auto id_value() -> ListenerId_t &;

    template<typename Callback>
    auto
    add_listener(const Callback &cb, ObserverMap<ListenerId_t, Callback> &callback_container) noexcept -> ListenerId;

    template<typename Callback>
    auto remove_listener(ListenerId_t id, ObserverMap<ListenerId_t, Callback> &callback_container) noexcept -> bool;

    template<typename Callback, typename ...Args>
    void notify(const ObserverMap<ListenerId_t, Callback> &callback_container, Args... args) const noexcept;

    template<typename KeyType, typename Callback>
    auto add_keyed_listener(KeyType key, const Callback &cb,
                            ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container) noexcept -> ListenerId;

    template<typename KeyType, typename Callback>
    auto remove_keyed_listener(KeyType key, ListenerId_t id,
                               ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container) noexcept -> bool;

    template<typename KeyType, typename Callback, typename ...Args>
    void keyed_notify(KeyType key, const ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container,
                      Args... args) const noexcept;

};

template<typename T>
auto Observable<T>::id_value() -> ListenerId_t & {
    static ListenerId_t the_id;
    return the_id;
}

template<typename T>
template<typename Callback>
auto Observable<T>::add_listener(const Callback &cb,
                                 ObserverMap<ListenerId_t, Callback> &callback_container) noexcept -> ListenerId {
    const auto id = ListenerId(++id_value());
    callback_container.emplace(id, cb);
    return id;
}

template<typename T>
template<typename Callback>
auto Observable<T>::remove_listener(const ListenerId_t id,
                                    ObserverMap<ListenerId_t, Callback> &callback_container) noexcept -> bool {
//    const auto it = callback_container.find(id);
//    if (it == callback_container.end()) {
//        return false;
//    }
    auto elements_removed = callback_container.erase(id);
    assert(elements_removed == 1);
    return elements_removed > 0;
}

template<typename T>
template<typename Callback, typename ...Args>
void Observable<T>::notify(const ObserverMap<ListenerId_t, Callback> &callback_container, Args... args) const noexcept {
    for (const auto &pair : callback_container) {
        pair.second(args...);
    }
}

template<typename T>
template<typename KeyType, typename Callback>
auto Observable<T>::add_keyed_listener(KeyType key, const Callback &cb,
                                       ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container) noexcept -> ListenerId {
    const auto id = ListenerId(++id_value());
    callback_container.at(key).emplace(id, cb);
    return id;
}

template<typename T>
template<typename KeyType, typename Callback>
auto Observable<T>::remove_keyed_listener(KeyType key, const ListenerId_t id,
                                          ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container) noexcept -> bool {
    auto cbs = callback_container.at(key);
//    const auto it = cbs.find(id);
//    if (it == cbs.end()) {
//        return false;
//    }
    auto elements_removed = cbs.erase(id);
    assert(elements_removed == 1);
    return elements_removed > 0;
}

template<typename T>
template<typename KeyType, typename Callback, typename ...Args>
void Observable<T>::keyed_notify(KeyType key,
                                 const ObserverMap<KeyType, ObserverMap<ListenerId_t, Callback>> &callback_container,
                                 Args... args) const noexcept {
    auto cbs = callback_container.at(key);
    for (const auto &pair : cbs) {
        pair.second(args...);
    }
}

#endif //TRANSMISSION_NETWORKS_APP_OBSERVABLE_H
