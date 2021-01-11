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

    template<typename ElementType>
    class OrderDerivedParentSet : public Computation<containers::ParentSet<ElementType>>,
                                  public abstract::Observable<OrderDerivedParentSet<ElementType>>,
                                  public abstract::Cacheable<OrderDerivedParentSet<ElementType>>,
                                  public abstract::Checkpointable<OrderDerivedParentSet<ElementType>, containers::ParentSet<ElementType>> {

        using ElementAddedCallback = std::function<void(ElementType *element)>;
        using ElementRemovedCallback = std::function<void(ElementType *element)>;

    protected:
        CREATE_EVENT(element_added, ElementAddedCallback)
        CREATE_EVENT(element_removed, ElementRemovedCallback)

    public:
        explicit OrderDerivedParentSet(parameters::Ordering<ElementType> &ordering, ElementType &child);
        containers::ParentSet<ElementType> value() noexcept override;
        void addDisallowedParent(ElementType* p);
        void addDisallowedParents(std::vector<ElementType*> p);
        void removeDisallowedParent(ElementType* p);
        void serialize() noexcept;

    protected:
        friend class abstract::Checkpointable<OrderDerivedParentSet<ElementType>, containers::ParentSet<ElementType>>;
        friend class abstract::Cacheable<OrderDerivedParentSet<ElementType>>;

        std::set<ElementType *> disallowedParents_{};
        parameters::Ordering<ElementType> &ordering_;
        ElementType &child_;
    };


    template<typename ElementType>
    OrderDerivedParentSet<ElementType>::OrderDerivedParentSet(parameters::Ordering<ElementType> &ordering, ElementType &child) : ordering_(ordering), child_(child) {
        ordering_.registerCacheableCheckpointTarget(this);

        ordering_.add_keyed_moved_left_listener(&child_, [=, this](ElementType *element) {
            if (!disallowedParents_.contains(element)) {
                this->setDirty();
                this->value_.insert(element);
                notify_element_added(element);
            }
        });

        ordering_.add_keyed_moved_right_listener(&child_, [=, this](ElementType *element) {
            if (!disallowedParents_.contains(element)) {
                this->setDirty();
                this->value_.erase(element);
                notify_element_removed(element);
            }
        });

        // Initialize the current parent set from the ordering
        for (auto &el : ordering.value()) {
            if (el != &child_) {
                this->value_.insert(el);
            } else {
                this->setClean();
                break;
            }
        }

        this->setDirty();
    }

    template<typename ElementType>
    void computation::OrderDerivedParentSet<ElementType>::addDisallowedParent(ElementType *p) {
        disallowedParents_.insert(p);
        if(this->value_.contains(p)) {
            this->value_.erase(p);
        }
    }

    template<typename ElementType>
    void computation::OrderDerivedParentSet<ElementType>::addDisallowedParents(std::vector<ElementType*> p) {
        for (const auto el : p) {
           addDisallowedParent(el);
        }
    }

    template<typename ElementType>
    void computation::OrderDerivedParentSet<ElementType>::removeDisallowedParent(ElementType *p) {
        disallowedParents_.erase(p);
    }

    template<typename ElementType>
    containers::ParentSet<ElementType> computation::OrderDerivedParentSet<ElementType>::value() noexcept {
        this->setClean();
        return this->value_;
    }

    template<typename ElementType>
    void computation::OrderDerivedParentSet<ElementType>::serialize() noexcept {
        std::cout << "{ ";
        for (auto &p : this->value()) {
            std::cout << *p << ", ";
        }
        std::cout << "}" << std::endl;
    }
}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
