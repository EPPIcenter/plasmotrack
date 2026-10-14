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

#include <algorithm>
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
        // The parent set itself rather than a copy, valid until it next changes.
        const containers::ParentSet<ElementType>& members() noexcept {
            this->setClean();
            return this->value_;
        }
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
        // Told of each allowed parent's changes: a parent's genotype change reaches only the parent sets that could
        // contain it, rather than every parent set through the ordering. A genotype change leaves the members as they
        // are, so this parent set is not saved with it; what depends on the members' genotypes saves with them itself.
        if (p != child_ and allowedParents_.insert(p).second) {
            if constexpr (requires { p->add_post_change_listener([]() {}); }) {
                // By address: a shared_ptr held in the parent's own listener would keep it alive.
                p->add_post_change_listener([this, parent = p.get()]() {
                    const auto found = std::ranges::lower_bound(this->value_, parent, containers::UidLess{}, [](const auto& element) { return element.get(); });
                    if (found != this->value_.end() and found->get() == parent) {
                        this->setDirty();
                        this->notify_element_changed(*found);
                    }
                });
            }
        }
    }

    template<typename ElementType, typename OrderingImpl>
    void OrderDerivedParentSet<ElementType, OrderingImpl>::addAllowedParents(const std::vector<std::shared_ptr<ElementType>>& p) {
        for (auto el : p) {
            addAllowedParent(el);
        }
    }


    template<typename ElementType, typename OrderingImpl>
    containers::ParentSet<ElementType> OrderDerivedParentSet<ElementType, OrderingImpl>::value() noexcept {
        this->setClean();
        return this->value_;
    }

}// namespace transmission_nets::core::computation


#endif//TRANSMISSION_NETWORKS_APP_ORDERDERIVEDPARENTSET_H
