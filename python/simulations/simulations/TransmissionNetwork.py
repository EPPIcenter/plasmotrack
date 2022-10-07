import numpy as np
from collections import defaultdict, deque
from simulations.Node import NodeLike, SourcePopulation
from typing import Callable, List, Optional


class TransmissionNetwork:
    def __init__(
            self,
            r0: float,
            node_constructor: Callable[[Optional[NodeLike], Optional[str]], NodeLike],
            offspring_sampler: Callable[..., int],
            source_population: SourcePopulation,
            num_founders: int,
            false_positive_rate=0.01,
            false_negative_rate=0.05,
            infection_duration_shape=10,
            infection_duration_scale=10,
            loss_rate=0.1,
            mutation_rate=0.0,
    ):
        self.num_founders = 0
        self.num_nodes = 0
        self.founders = []
        self.nodes = []
        self.edges = defaultdict(list)
        self.edge_list = []
        self.parent_list = deque()

        self.r0 = r0
        self.node_constructor = node_constructor
        self.offspring_sampler = offspring_sampler
        self.source_population = source_population
        self.num_founders = num_founders
        self.infection_duration_shape = infection_duration_shape
        self.infection_duration_scale = infection_duration_scale
        self.false_positive_rate = false_positive_rate
        self.false_negative_rate = false_negative_rate
        self.loss_rate = loss_rate
        self.mutation_rate = mutation_rate

        # construct the founders within the network
        for _ in range(num_founders):
            self.add_founder(self.source_population)

        # generate transmissions at rate r0
        while self.parent_list:
            if len(self.parent_list) > 1:
                num_parents = 1 if np.random.random(1) > .2 else 2
            else:
                num_parents = 1
            parents = [self.parent_list.popleft() for _ in range(num_parents)]
            self.add_children(parents)

    def add_edge(self, source: NodeLike, dest: NodeLike):
        self.edges[source].append(dest)

    def add_founder(self, source: SourcePopulation):
        node = self.node_constructor(
            false_positive_rate=self.false_positive_rate,
            false_negative_rate=self.false_negative_rate,
            infection_duration_scale=self.infection_duration_scale,
            infection_duration_shape=self.infection_duration_shape,
            mutation_rate=self.mutation_rate,
            loss_rate=self.loss_rate,
            source=source,
        )
        self.founders.append(node)
        self.nodes.append(node)
        self.edges[node] = []
        self.parent_list.append(node)
        self.num_founders += 1
        self.num_nodes += 1

    def add_children(self, parents):
        num_offspring = self.offspring_sampler(self.r0)
        for _ in range(num_offspring):
            child = self.node_constructor(
                false_positive_rate=self.false_positive_rate,
                false_negative_rate=self.false_negative_rate,
                infection_duration_scale=self.infection_duration_scale,
                infection_duration_shape=self.infection_duration_shape,
                mutation_rate=self.mutation_rate,
                loss_rate=self.loss_rate,
                parents=parents
            )
            self.nodes.append(child)
            self.edges[child] = []
            for parent in parents:
                self.add_edge(parent, child)
                self.edge_list.append((parent, child))
            self.parent_list.append(child)
            self.num_nodes += 1


def subsample_network(tnet: TransmissionNetwork, prop: float):
    rng = np.random.default_rng()
    nodes_subset = []
    nodes_to_remove = []
    edges = {
        k: list(tnet.edges[k]) for k in tnet.edges.keys()
    }  # create copy which we will modify by removing edges
    for node in tnet.nodes:
        if rng.binomial(1, prop):
            nodes_subset.append(node)
        else:
            nodes_to_remove.append(node)
            child_set = edges.pop(node)
            parent_set = [
                parent_node
                for parent_node in edges.keys()
                if node in edges[parent_node]
            ]
            for parent_node in parent_set:
                for child_node in child_set:
                    edges[parent_node].append(child_node)
    return nodes_subset, edges


def reduce_network(tnet: TransmissionNetwork, nodes: List[NodeLike]):
    nodes_subset = []
    nodes_to_remove = []
    edges = {
        k: list(tnet.edges[k]) for k in tnet.edges.keys()
    }  # create copy which we will modify by removing edges
    for node in tnet.nodes:
        if node in nodes:
            nodes_subset.append(node)
        else:
            nodes_to_remove.append(node)
            child_set = edges.pop(node)
            parent_set = [
                parent_node
                for parent_node in edges.keys()
                if node in edges[parent_node]
            ]
            if (
                    parent_set
            ):  # does the node have a parent? if so, make that the parent of the children nodes
                for parent_node in parent_set:
                    for child_node in child_set:
                        edges[parent_node].append(child_node)
            else:  # otherwise, deeply connect all the children nodes
                for i in range(len(child_set)):
                    for j in range(i + 1, len(child_set)):
                        edges[child_set[i]].append(child_set[j])
                        edges[child_set[j]].append(child_set[i])

    return nodes_subset, edges
