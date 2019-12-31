//
// Created by Maxwell Murphy on 12/8/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_NODE_H
#define TRANSMISSION_NETWORKS_APP_NODE_H

#include <iostream>
#include <vector>
#include <optional>
#include <variant>
#include "AbstractNode.h"

template<typename T>
class Node : public AbstractNode {
protected:
    bool is_dirty_{true};

    T value_;
    std::optional<T> saved_state_;

    std::vector<AbstractNode *> child_dependencies_{}; // Nodes that depend on this instance
    boost::container::flat_set<AbstractNode *> dirty_dependencies_{}; // Dirty Nodes this Node depends on

public:

    virtual T value() = 0;

    virtual T peek() = 0;

    bool isDirty();

    void setDirty();

    void setDirty(AbstractNode* ptr);

    void setClean();

    void saveState();

    void restoreState();

    void acceptState();

    void addDependency(AbstractNode* ptr);
};


template<typename T>
void Node<T>::setDirty() { // Propagate dirty through the dependency tree
    if (is_dirty_) {
        return;
    }

    is_dirty_ = true;
    for (auto node : child_dependencies_) {
        node->setDirty(this);
    }
}

template<typename T>
void Node<T>::setDirty(AbstractNode* ptr) {
    dirty_dependencies_.insert(ptr);

    if (is_dirty_) {
        return;
    }

    is_dirty_ = true;
    for (auto node : child_dependencies_) {
        node->setDirty(this);
    }
}

template<typename T>
void Node<T>::setClean() {
    is_dirty_ = false;
    dirty_dependencies_.clear();
}

template<typename T>
void Node<T>::saveState() {
    assert(!is_dirty_); // graph must be clean to save state.

    if(saved_state_) {
        // early break because already saved.
        return;
    }

    for (auto node : child_dependencies_) {
        node->saveState();
    }
    saved_state_ = T(value_);
}

template<typename T>
void Node<T>::restoreState() {
    if (!saved_state_) {
        // early break because already restored.
        return;
    }

    value_ = *saved_state_;
    this->setClean();
    for (auto node : child_dependencies_) {
        node->restoreState();
    }
    saved_state_.reset();
}

template<typename T>
void Node<T>::acceptState() {
    if (!saved_state_) {
        // early break because already accepted.
        return;
    }

    this->setClean();
    for (auto node : child_dependencies_) {
        node->acceptState();
    }
    saved_state_.reset();
}

template<typename T>
bool Node<T>::isDirty() {
    return is_dirty_;
}

template<typename T>
void Node<T>::addDependency(AbstractNode* ptr) {
    child_dependencies_.push_back(ptr);
}

#endif //TRANSMISSION_NETWORKS_APP_NODE_H
