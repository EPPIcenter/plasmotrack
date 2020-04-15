//
// Created by Maxwell Murphy on 2/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERING_H
#define TRANSMISSION_NETWORKS_APP_ORDERING_H

#include <iostream>
#include <vector>

#include "core/parameters/Parameter.h"

template <typename T>
class Ordering : public Parameter<std::vector<T*>> {

    using MovedCallback = std::function<void(T* element)>;
    CREATE_KEYED_EVENT(moved_left, T*, MovedCallback) // Notifies that an element has been moved left of key
    CREATE_KEYED_EVENT(moved_right, T*, MovedCallback) // Notifies that an element has been moved right of key

public:
    explicit Ordering() noexcept;

    explicit Ordering(std::vector<T*> refs) noexcept;

    void swap(int a, int b) noexcept;

    void addElement(T* ref) noexcept;

    void addElements(std::vector<T*> refs) noexcept;

    friend std::ostream &operator<<(std::ostream &os, const Ordering &list) noexcept {
        for (unsigned long i = 0; i < list.value_.size(); ++i) {
            os << "Element " << i << ": " << *list.value_[i] << std::endl;
        }
        return os;
    };

private:
    void notifySwap(int left_idx, int right_idx) noexcept;
};

template<typename T>
Ordering<T>::Ordering() noexcept : Parameter<std::vector<T*>>() {}

template<typename T>
Ordering<T>::Ordering(std::vector<T *> refs) noexcept {
    addElements(refs);
}

template<typename T>
void Ordering<T>::swap(int a, int b) noexcept {
    if(a != b) {
        T* tmp = this->value_.at(a);
        this->value_.at(a) = this->value_.at(b);
        this->value_.at(b) = tmp;

        if (a < b) {
            notifySwap(a, b);
        } else {
            notifySwap(b, a);
        }
    }
}

template<typename T>
void Ordering<T>::addElement(T *ref) noexcept {
    this->value_.push_back(ref);
    register_moved_left_listener_key(ref);
    register_moved_right_listener_key(ref);
}

template<typename T>
void Ordering<T>::addElements(std::vector<T *> refs) noexcept {
    for (auto& ref: refs) {
        addElement(ref);
    }
}

template<typename T>
void Ordering<T>::notifySwap(const int left_idx, const int right_idx) noexcept {
    keyed_notify_moved_right(this->value_[left_idx], this->value_[right_idx]);
    keyed_notify_moved_left(this->value_[right_idx], this->value_[left_idx]);
    for (int i = left_idx + 1; i < right_idx; ++i) {
        keyed_notify_moved_left(this->value_[right_idx], this->value_[i]);
        keyed_notify_moved_right(this->value_[left_idx], this->value_[i]);
        keyed_notify_moved_right(this->value_[i], this->value_[right_idx]);
        keyed_notify_moved_left(this->value_[i], this->value_[left_idx]);
    }
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERING_H
