//
// Created by Maxwell Murphy on 10/7/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDORDERDERIVEDPARENTSET_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDORDERDERIVEDPARENTSET_H

#include <set>

#include "core/computation/OrderDerivedParentSet.h"

namespace transmission_nets::core::computation {

    template<typename ElementType>
    class ConstrainedOrderDerivedParentSet : public OrderDerivedParentSet<ElementType> {

    public:
        ConstrainedOrderDerivedParentSet(parameters::Ordering<ElementType>& ordering, std::set<ElementType*> allowedParents, ElementType& child);

    private:
        std::set<ElementType*> allowedParents_;
    };

    template<typename ElementType>
    ConstrainedOrderDerivedParentSet<ElementType>::ConstrainedOrderDerivedParentSet(parameters::Ordering<ElementType>& ordering, std::set<ElementType*> allowedParents, ElementType& child) : ordering_(ordering), allowedParents_(allowedParents), child_(child) {
        this->ordering_.registerCacheableCheckpointTarget(this);

        this->ordering_.add_keyed_moved_left_listener(&(this->child_), [=, this](ElementType* element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.insert(element);
                this->notify_element_added(element);
            }
        });

        this->ordering_.add_keyed_moved_right_listener(&(this->child_), [=, this](ElementType* element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.erase(element);
                this->notify_element_removed(element);
            }
        });

        for (auto& el : ordering.value()) {
            if (el != &(this->child_)) {
                if (allowedParents.contains(el)) {
                    this->value_.insert(el);
                }
            }
        }
    }


}// namespace transmission_nets::core::computation

#endif//TRANSMISSION_NETWORKS_APP_CONSTRAINEDORDERDERIVEDPARENTSET_H
