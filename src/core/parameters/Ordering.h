//
// Created by Maxwell Murphy on 2/1/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERING_H
#define TRANSMISSION_NETWORKS_APP_ORDERING_H


#include "core/parameters/Parameter.h"

#include <iostream>
#include <memory>
#include <vector>


namespace transmission_nets::core::parameters {
    template<typename T>
    class Ordering : public Parameter<std::vector<std::shared_ptr<T>>> {

        using MovedCallback = std::function<void(std::shared_ptr<T> element)>;
        using ChangedCallback = std::function<void(std::shared_ptr<T> element)>;
        CREATE_KEYED_EVENT(moved_left, std::shared_ptr<T>, MovedCallback) // Notifies that an element has been moved left of key
        CREATE_KEYED_EVENT(moved_right, std::shared_ptr<T>, MovedCallback)// Notifies that an element has been moved right of key
        CREATE_EVENT(element_changed, ChangedCallback);


    public:
        explicit Ordering() noexcept;

        explicit Ordering(std::vector<std::shared_ptr<T>> refs) noexcept;

        void swap(int a, int b) noexcept;

        void addElement(std::shared_ptr<T> ref) noexcept;

        void addElements(const std::vector<std::shared_ptr<T>>& refs) noexcept;

        friend std::ostream& operator<<(std::ostream& os, const Ordering& list) noexcept {
            for (unsigned long i = 0; i < list.value_.size(); ++i) {
                os << "Element " << i << ": " << *list.value_[i] << "\n";
            }
            return os;
        };


    private:
        void notifySwap(int left_idx, int right_idx) noexcept;
        void elementChanged(std::shared_ptr<T> ref) noexcept;
    };

    template<typename T>
    Ordering<T>::Ordering() noexcept : Parameter<std::vector<std::shared_ptr<T>>>() {}

    template<typename T>
    Ordering<T>::Ordering(std::vector<std::shared_ptr<T>> refs) noexcept {
        addElements(refs);
    }

    template<typename T>
    void Ordering<T>::swap(int a, int b) noexcept {
        if (a != b) {
            this->notify_pre_change();
            auto tmp           = this->value_.at(a);
            this->value_.at(a) = this->value_.at(b);
            this->value_.at(b) = tmp;

            if (a < b) {
                notifySwap(a, b);
            } else {
                notifySwap(b, a);
            }
            this->notify_post_change();
        }
    }

    template<typename T>
    void Ordering<T>::addElement(std::shared_ptr<T> ref) noexcept {
        this->notify_pre_change();
        this->value_.push_back(ref);
        register_moved_left_listener_key(ref);
        register_moved_right_listener_key(ref);
        this->notify_post_change();
    }

    template<typename T>
    void Ordering<T>::addElements(const std::vector<std::shared_ptr<T>>& refs) noexcept {
        for (auto& ref : refs) {
            addElement(std::move(ref));
        }
    }

    template<typename T>
    void Ordering<T>::notifySwap(const int left_idx, const int right_idx) noexcept {
        keyed_notify_moved_right(this->value_[left_idx], this->value_[right_idx]);
        keyed_notify_moved_left(this->value_[right_idx], this->value_[left_idx]);
        for (int i = left_idx + 1; i < right_idx; ++i) {

            // value_[i] has been moved left of element at value_[right_idx]
            keyed_notify_moved_left(this->value_[right_idx], this->value_[i]);

            // value_[right_idx] has been moved right of element at value_[i]
            keyed_notify_moved_right(this->value_[i], this->value_[right_idx]);

            // value_[i] has been moved right of element at value_[left_idx]
            keyed_notify_moved_right(this->value_[left_idx], this->value_[i]);

            // value_[left_idx] has been moved left of element at value_[i]
            keyed_notify_moved_left(this->value_[i], this->value_[left_idx]);
        }
    }

    template<typename T>
    void Ordering<T>::elementChanged(std::shared_ptr<T> ref) noexcept {
        this->notify_element_changed(ref);
        this->setDirty();
    }

}// namespace transmission_nets::core::parameters


#endif//TRANSMISSION_NETWORKS_APP_ORDERING_H
