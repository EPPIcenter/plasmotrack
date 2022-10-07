from collections import defaultdict
from typing import Dict, List

"""
This is a network topology implementation that can be utilized to emulate various epidemiological settings. Transmission events within a population can be represented
either as a DAG (where the node is an infection event) or as a digraph (where the node is the individual). The digraph can be thought of as a flattening of the DAG along the temporal
axis.
"""


class Node:
    _nodes: Dict[str, "Node"] = {}

    def __init__(self, label: str, epi_data=None) -> None:
        self.label = label
        self.epi_data = epi_data or {}
        Node._nodes[self.label] = self

    @classmethod
    def get_node(cls, label: str) -> "Node":
        return cls._nodes[label]

    @classmethod
    def get_node_labels(cls) -> List[str]:
        return cls._nodes.keys()


class Network:
    def __init__(self) -> None:
        self.child_sets: Dict[str, List[str]] = defaultdict(list)
        self.parent_sets: Dict[str, List[str]] = defaultdict(list)
        self.total_nodes = 0

    def add_edge(self, source: str, dest: str):
        if source not in self.child_sets and source not in self.parent_sets:
            self.total_nodes += 1
        if dest not in self.parent_sets and dest not in self.child_sets:
            self.total_nodes += 1
        self.child_sets[source].append(dest)
        self.parent_sets[dest].append(source)

    def get_children(self, source: str):
        return self.child_sets[source]

    def get_parents(self, source: str):
        return self.parent_sets[source]

    def get_roots(self):
        return [
            node_label
            for node_label in Node.get_node_labels()
            if not self.parent_sets[node_label]
        ]

    def get_topological_sorting(self):
        """
        Returns a topological ordering of the nodes in the network. This is useful when needing to process the nodes
        in an ordering such that all parents are processed before their children.
        """
        t = []
        all_nodes = sorted(
            list(
                set(
                    list(self.child_sets.keys())
                    + self.get_roots()
                    + list(self.parent_sets.keys())
                )
            ),
            key=lambda x: int(x),
        )
        total_nodes = len(all_nodes)
        in_degree = [0] * total_nodes
        visited = [False] * total_nodes

        for i in range(total_nodes):
            in_degree[i] = len(self.get_parents(all_nodes[i]))

        print("In Degree: ", in_degree)

        to_process = []
        for i in range(total_nodes):
            if in_degree[i] == 0:
                to_process.append(i)
                visited[i] = True

        while to_process:
            tar_index = to_process.pop(0)
            # add the node to the topological ordering
            t.append(all_nodes[tar_index])

            # for each child of the node, decrement the in-degree
            for j in range(total_nodes):
                if (
                    all_nodes[j] in self.get_children(all_nodes[tar_index])
                    and visited[j] is False
                ):
                    in_degree[j] -= 1
                    if in_degree[j] == 0:
                        to_process.append(j)
                        visited[j] = True

        return t
