//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
#define TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H

#include "core/abstract/observables/Cacheable.h"
#include "core/abstract/observables/Checkpointable.h"
#include "core/abstract/observables/Observable.h"
#include "core/computation/Computation.h"
#include "core/containers/ParentSet.h"
#include "core/parameters/Ordering.h"

#include <memory>
#include <set>

namespace transmission_nets::core::computation {

    template<typename ElementType, typename OrderingImpl>
    class OrderDerivedParentSet : public Computation<containers::ParentSet<ElementType>>,
                                  public abstract::Observable<OrderDerivedParentSet<ElementType, OrderingImpl>>,
                                  public abstract::Cacheable<OrderDerivedParentSet<ElementType, OrderingImpl>>,
                                  public abstract::Checkpointable<OrderDerivedParentSet<ElementType, OrderingImpl>, containers::ParentSet<ElementType>> {

        using ElementAddedCallback   = std::function<void(std::shared_ptr<ElementType> element)>;
        using ElementRemovedCallback = std::function<void(std::shared_ptr<ElementType> element)>;
        using ElementChangedCallback = std::function<void(std::shared_ptr<ElementType> element)>;

    protected:
        CREATE_EVENT(element_added, ElementAddedCallback)
        CREATE_EVENT(element_removed, ElementRemovedCallback)
        CREATE_EVENT(element_changed, ElementChangedCallback)

    public:
        explicit OrderDerivedParentSet(std::shared_ptr<OrderingImpl> ordering,
                                       std::shared_ptr<ElementType> child,
                                       const std::vector<std::shared_ptr<ElementType>>& allowedParents = {});
        containers::ParentSet<ElementType> value() noexcept override;
        void addAllowedParent(std::shared_ptr<ElementType> p);
        void addAllowedParents(const std::vector<std::shared_ptr<ElementType>>& p);

        std::shared_ptr<OrderingImpl> ordering_;
        std::shared_ptr<ElementType> child_;
        std::set<std::shared_ptr<ElementType>> allowedParents_{};

    protected:
        friend class abstract::Checkpointable<OrderDerivedParentSet, containers::ParentSet<ElementType>>;
        friend class abstract::Cacheable<OrderDerivedParentSet>;

    };

    template<typename ElementType, typename OrderingImpl>
    OrderDerivedParentSet<ElementType, OrderingImpl>::OrderDerivedParentSet(std::shared_ptr<OrderingImpl> ordering, std::shared_ptr<ElementType> child, const std::vector<std::shared_ptr<ElementType>>& allowedParents) : ordering_(std::move(ordering)), child_(std::move(child)) {
        ordering_->registerCacheableCheckpointTarget(this);

        // If the element moves to the left of child_, trigger this event
        ordering_->add_keyed_moved_left_listener(child_, [=, this](std::shared_ptr<ElementType> element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.insert(element);
                this->notify_element_added(element);
            }
        });

        // If the element moves to the right of child_, trigger this event
        ordering_->add_keyed_moved_right_listener(child_, [=, this](std::shared_ptr<ElementType> element) {
            if (allowedParents_.contains(element)) {
                this->setDirty();
                this->value_.erase(element);
                this->notify_element_removed(element);
            }
        });

        // If an element changes, trigger this event
        ordering_->add_element_changed_listener([=, this](std::shared_ptr<ElementType> element) {
            if (this->value_.contains(element)) {
                this->setDirty();
                this->notify_element_changed(element);
            }
        });

        this->addAllowedParents(allowedParents);

        // Initialize the current parent set from the ordering
        for (auto& el : ordering_->value()) {
            if (el != child_) {
                if (allowedParents_.contains(el)) {
                    this->value_.insert(el);
                }
            } else {
                break;
            }
        }
        this->setDirty();
    }

    template<typename ElementType, typename OrderingImpl>
    void OrderDerivedParentSet<ElementType, OrderingImpl>::addAllowedParent(std::shared_ptr<ElementType> p) {
        if (p != child_) {
            allowedParents_.insert(p);
        }
    }

    template<typename ElementType, typename OrderingImpl>
    void OrderDerivedParentSet<ElementType, OrderingImpl>::addAllowedParents(const std::vector<std::shared_ptr<ElementType>>& p) {
        for (auto el : p) {
            addAllowedParent(std::move(el));
        }
    }


    template<typename ElementType, typename OrderingImpl>
    containers::ParentSet<ElementType> OrderDerivedParentSet<ElementType, OrderingImpl>::value() noexcept {
        this->setClean();
        return this->value_;
    }

}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
