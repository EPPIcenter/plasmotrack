import json
from collections import defaultdict, deque
from typing import Callable, Optional

import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np

from simulations.Node import NodeLike, SimpleNode, SourcePopulation


class TransmissionNetwork:
    def __init__(
        self,
        r0: float,
        node_constructor: Callable[[Optional[NodeLike], Optional[str]], NodeLike],
        offspring_sampler: Callable[..., int],
        false_positive_rate=0.01,
        false_negative_rate=0.1,
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

        self.node_constructor = node_constructor
        self.offspring_sampler = offspring_sampler
        self.r0 = r0
        self.false_positive_rate = false_positive_rate
        self.false_negative_rate = false_negative_rate
        self.loss_rate = loss_rate
        self.mutation_rate = mutation_rate

    def add_edge(self, source: NodeLike, dest: NodeLike):
        self.edges[source].append(dest)

    def add_founder(self, source: SourcePopulation):
        node = self.node_constructor(
            false_positive_rate=self.false_positive_rate,
            false_negative_rate=self.false_negative_rate,
            mutation_rate=self.mutation_rate,
            loss_rate=self.loss_rate,
            source=source,
        )
        self.founders.append(node)
        self.nodes.append(node)
        self.edges[node] = []
        # self.edge_list.append((source, node))
        self.parent_list.append(node)
        self.num_founders += 1
        self.num_nodes += 1

    def add_children(self, parent):
        num_offspring = self.offspring_sampler(self.r0)
        for _ in range(num_offspring):
            child = self.node_constructor(
                false_positive_rate=self.false_positive_rate,
                false_negative_rate=self.false_negative_rate,
                mutation_rate=self.mutation_rate,
                loss_rate=self.loss_rate,
                parent=parent,
            )
            self.nodes.append(child)
            self.edges[child] = []
            self.add_edge(parent, child)
            self.edge_list.append((parent, child))
            self.parent_list.append(child)
            self.num_nodes += 1


def calculate_distance(node1, node2) -> int:
    d = 0
    for locus in node1.observed_alleles.keys():
        node1_a = node1.observed_alleles[locus]
        node2_a = node2.observed_alleles[locus]
        d += np.count_nonzero(node1_a != node2_a)
    return d


rng = np.random.default_rng()
t1 = TransmissionNetwork(
    r0=0.8, node_constructor=SimpleNode, offspring_sampler=lambda x: rng.poisson(x)
)

allele_frequencies = {}
for k in range(20):
    allele_frequencies[f"L{k+1}"] = rng.dirichlet(np.array([1] * 10))

source1 = SourcePopulation(allele_frequencies=allele_frequencies, mean_coi=4, label="S")


for i in range(25):
    t1.add_founder(source1)

while t1.parent_list:
    parent = t1.parent_list.popleft()
    t1.add_children(parent)


pairwise_distances = defaultdict(list)
dist_list = []

for i in range(len(t1.nodes)):
    for j in range(i + 1, len(t1.nodes)):
        node_i = t1.nodes[i]
        node_j = t1.nodes[j]
        d = calculate_distance(node_i, node_j)
        pairwise_distances[node_i].append((node_j, d))
        pairwise_distances[node_j].append((node_i, d))
        dist_list.append(d)


disallowed_parents = dict()
allowed_parents = dict()

for node in t1.nodes:
    dists = pairwise_distances[node]
    dists.sort(key=lambda x: x[1])
    disallowed_parents[node] = dists[25:]
    allowed_parents[node] = dists[:25]


in_allowed_parents = dict()
for node, ps in allowed_parents.items():
    parents = [_[0] for _ in ps]
    print(len(ps))
    print(len(parents))
    in_allowed_parents[node] = node.parent in parents


out = {
    "loci": [
        {
            "locus": locus,
            "allele_freqs": list(allele_freqs),
            "num_alleles": len(allele_freqs),
        }
        for (locus, allele_freqs) in allele_frequencies.items()
    ],
    "nodes": [
        {
            "id": node.label,
            "latent_genotype": [
                {"locus": locus, "genotype": "".join(map(str, list(genotype)))}
                for (locus, genotype) in node.true_alleles.items()
            ],
            "observed_genotype": [
                {"locus": locus, "genotype": "".join(map(str, list(genotype)))}
                for (locus, genotype) in node.observed_alleles.items()
            ],
            "disallowed_parents": [p[0].label for p in disallowed_parents[node]],
        }
        for node in t1.nodes
    ],
    "network": [
        {"from": parent.label, "to": child.label} for (parent, child) in t1.edge_list
    ],
}

with open(
    "/Users/maxwellmurphy/Workspace/transmission_nets/test/resources/JSON/nodes5.json",
    "w",
) as f:
    json.dump(out, f, separators=(",", ":"))
