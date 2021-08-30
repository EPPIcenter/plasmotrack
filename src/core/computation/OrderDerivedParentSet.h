//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
#define TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H

#include <iostream>
#include <set>

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Ordering.h"


namespace transmission_nets::core::computation {

    template<typename ElementType, typename OrderingImpl>
    class OrderDerivedParentSet : public Computation<containers::ParentSet<ElementType>>,
                                  public abstract::Observable<OrderDerivedParentSet<ElementType, OrderingImpl>>,
                                  public abstract::Cacheable<OrderDerivedParentSet<ElementType, OrderingImpl>>,
                                  public abstract::Checkpointable<OrderDerivedParentSet<ElementType, OrderingImpl>, containers::ParentSet<ElementType>> {

        using ElementAddedCallback = std::function<void(ElementType *element)>;
        using ElementRemovedCallback = std::function<void(ElementType *element)>;

    protected:
        CREATE_EVENT(element_added, ElementAddedCallback)
        CREATE_EVENT(element_removed, ElementRemovedCallback)

    public:
        explicit OrderDerivedParentSet(OrderingImpl *ordering, ElementType *child);
        containers::ParentSet<ElementType> value() noexcept override;
        void addAllowedParent(ElementType* p);
        void addAllowedParents(std::vector<ElementType*> p);
        void removeAllowedParent(ElementType* p);
        void serialize() noexcept;

    protected:
        friend class abstract::Checkpointable<OrderDerivedParentSet<ElementType, OrderingImpl>, containers::ParentSet<ElementType>>;
        friend class abstract::Cacheable<OrderDerivedParentSet<ElementType, OrderingImpl>>;

        std::set<ElementType *> allowedParents_{};
        OrderingImpl *ordering_;
        ElementType *child_;
    };


    template<typename ElementType, typename OrderingImpl>
    OrderDerivedParentSet<ElementType, OrderingImpl>::OrderDerivedParentSet(OrderingImpl *ordering, ElementType *child) : ordering_(ordering), child_(child) {
        ordering_->registerCacheableCheckpointTarget(this);

        ordering_->add_keyed_moved_left_listener(child_, [=, this](ElementType *element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.insert(element);
                notify_element_added(element);
            }
        });

        ordering_->add_keyed_moved_right_listener(child_, [=, this](ElementType *element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.erase(element);
                notify_element_removed(element);
            }
        });

        // Initialize the current parent set from the ordering
        for (auto &el : ordering_->value()) {
            if (el != child_ and allowedParents_.contains(el)) {
                this->value_.insert(el);
            } else {
                this->setClean();
                break;
            }
        }

        this->setDirty();
    }

    template<typename ElementType, typename OrderingImpl>
    void computation::OrderDerivedParentSet<ElementType, OrderingImpl>::addAllowedParent(ElementType *p) {
        allowedParents_.insert(p);
//        if(!(this->value_.contains(p))) {
//            this->value_.insert(p);
//        }
    }

    template<typename ElementType, typename OrderingImpl>
    void computation::OrderDerivedParentSet<ElementType, OrderingImpl>::addAllowedParents(std::vector<ElementType*> p) {
        for (const auto el : p) {
           addAllowedParent(el);
        }
    }

    template<typename ElementType, typename OrderingImpl>
    void computation::OrderDerivedParentSet<ElementType, OrderingImpl>::removeAllowedParent(ElementType *p) {
        allowedParents_.erase(p);
    }

    template<typename ElementType, typename OrderingImpl>
    containers::ParentSet<ElementType> computation::OrderDerivedParentSet<ElementType, OrderingImpl>::value() noexcept {
        this->setClean();
        return this->value_;
    }

    template<typename ElementType, typename OrderingImpl>
    void computation::OrderDerivedParentSet<ElementType, OrderingImpl>::serialize() noexcept {
        std::cout << "{ ";
        for (auto &p : this->value()) {
            std::cout << *p << ", ";
        }
        std::cout << "}" << std::endl;
    }
}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
