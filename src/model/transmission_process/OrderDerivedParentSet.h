//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
#define TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H

#include <iostream>

#include "core/computation/Computation.h"

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"

#include "core/parameters/Ordering.h"

#include "core/containers/ParentSet.h"

template<typename ElementType>
class OrderDerivedParentSet : public Computation<ParentSet<ElementType>>,
                              public Observable<OrderDerivedParentSet<ElementType>>,
                              public Cacheable<OrderDerivedParentSet<ElementType>>,
                              public Checkpointable<OrderDerivedParentSet<ElementType>, ParentSet<ElementType>> {

    using ElementAddedCallback = std::function<void(ElementType* element)>;
    using ElementRemovedCallback = std::function<void(ElementType* element)>;
    CREATE_EVENT(element_added, ElementAddedCallback);
    CREATE_EVENT(element_removed, ElementRemovedCallback);

public:
    explicit OrderDerivedParentSet(Ordering<ElementType> &ordering, ElementType &child) : ordering_(ordering), child_(child) {
        ordering_.registerCacheableCheckpointTarget(*this);

        ordering_.add_keyed_moved_left_listener(&child_, [&](ElementType* element) {
            this->value_.insert(element);
            this->notify_element_added(element);
            this->setDirty();
        });

        ordering_.add_keyed_moved_right_listener(&child_, [&](ElementType* element) {
            this->value_.erase(element);
            this->notify_element_removed(element);
            this->setDirty();
        });

        // Initialize the current parent set from the ordering
        for(auto &el : ordering.value()) {
            if (el != &child_) {
                this->value_.insert(el);
            } else {
                this->setClean();
                break;
            }
        }
    };

    ParentSet<ElementType> value() noexcept override {
        this->setClean();
        return this->value_;
    };

    void printSet() noexcept {
        std::cout << "{ ";
        for (auto &p : this->value()) {
            std::cout << *p << ", ";
        }
        std::cout << "}" << std::endl;
    }

private:
    friend class Checkpointable<OrderDerivedParentSet<ElementType>, ParentSet<ElementType>>;
    friend class Cacheable<OrderDerivedParentSet<ElementType>>;

    Ordering<ElementType> &ordering_;
    ElementType &child_;
};

#endif //TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
