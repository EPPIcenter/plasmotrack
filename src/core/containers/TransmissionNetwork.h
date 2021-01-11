//
// Created by Maxwell Murphy on 2/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GRAPH_H
#define TRANSMISSION_NETWORKS_APP_GRAPH_H

#include <boost/container/flat_map.hpp>

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"

#include "core/parameters/Parameter.h"

#include "core/containers/ParentSet.h"


template<typename NodeValueImpl>
class TransmissionNetwork :  public Observable<TransmissionNetwork<NodeValueImpl>>,
                             public UncacheablePassthrough<TransmissionNetwork<NodeValueImpl>>,
                             public CheckpointablePassthrough<TransmissionNetwork<NodeValueImpl>> {


public:
    TransmissionNetwork() = default;

    void addNode(NodeValueImpl *node) noexcept;

    void addNodes(std::vector<NodeValueImpl*> nodes) noexcept;

    void addEdge(NodeValueImpl* parent, NodeValueImpl* child) noexcept;

    void removeEdge(NodeValueImpl* parent, NodeValueImpl* child) noexcept;

    bool createsCycle(NodeValueImpl* parent, NodeValueImpl* child) noexcept;

    Parameter<ParentSet<NodeValueImpl>>* parentSet(NodeValueImpl* child) noexcept;

    [[nodiscard]] unsigned int totalNodes() const noexcept;

    const std::vector<NodeValueImpl*>& nodes() const noexcept;

    std::string serialize() const noexcept;

private:
    boost::container::flat_map<NodeValueImpl*, Parameter<ParentSet<NodeValueImpl>>*> parentSets_{};
    std::vector<NodeValueImpl*> nodes_{};
};

template<typename NodeValueImpl>
void TransmissionNetwork<NodeValueImpl>::addNode(NodeValueImpl *node) noexcept {
    parentSets_.emplace(node, new Parameter<ParentSet<NodeValueImpl>>());
    nodes_.push_back(node);
}

template<typename NodeValueImpl>
void TransmissionNetwork<NodeValueImpl>::addNodes(std::vector<NodeValueImpl*> nodes) noexcept {
    for (auto node : nodes) {
        addNode(node);
    }
}

template<typename NodeValueImpl>
void TransmissionNetwork<NodeValueImpl>::addEdge(NodeValueImpl* parent, NodeValueImpl* child) noexcept {
    assert(parentSet(child)->value().find(parent) == parentSet(child)->value().end());
    auto tmpParentSet = parentSet(child)->value();
    tmpParentSet.insert(parent);
    parentSet(child)->setValue(tmpParentSet);
}


template<typename NodeValueImpl>
void TransmissionNetwork<NodeValueImpl>::removeEdge(NodeValueImpl* parent, NodeValueImpl* child) noexcept {
    assert(parentSet(child)->value().find(parent) != parentSet(child)->value().end());
    auto tmpParentSet = parentSet(child)->value();
    tmpParentSet.erase(parent);
    parentSet(child)->setValue(tmpParentSet);
}

template<typename NodeValueImpl>
bool TransmissionNetwork<NodeValueImpl>::createsCycle(NodeValueImpl* parent,
                                                      NodeValueImpl* child) noexcept {
    // is potential child topologically before parent

    if(parent == child) {
        return true;
    }

    auto tmpParentTracker_ = parentSet(parent)->value();
    while(tmpParentTracker_.size() > 0) {
        auto el = tmpParentTracker_.begin();
        if (*el == child) {
            return true;
        } else {
            tmpParentTracker_.erase(el);
            tmpParentTracker_.insert(parentSet(*el)->value().begin(), parentSet(*el)->value().end());
        }
    }
    return false;
}


template<typename NodeValueImpl>
Parameter<ParentSet<NodeValueImpl>>* TransmissionNetwork<NodeValueImpl>::parentSet(NodeValueImpl* child) noexcept {
    return parentSets_.at(child);
}

template<typename NodeValueImpl>
unsigned int TransmissionNetwork<NodeValueImpl>::totalNodes() const noexcept {
    return parentSets_.size();
}

template<typename NodeValueImpl>
const std::vector<NodeValueImpl *>& TransmissionNetwork<NodeValueImpl>::nodes() const noexcept {
    return nodes_;
}

template<typename NodeValueImpl>
std::string TransmissionNetwork<NodeValueImpl>::serialize() const noexcept {
    std::string out;
    for (const auto [child, parent_set_param] : parentSets_) {
        const auto& parent_set = parent_set_param->value();
        if (parent_set.size() == 0) {
            out += "S-" + child->serialize() + ",";
        } else {
            for (const auto parent : parent_set) {
                out += parent->serialize() + "-" + child->serialize() + ",";
            }
        }
    }
    out.pop_back();
    return out;
}


#endif //TRANSMISSION_NETWORKS_APP_GRAPH_H
