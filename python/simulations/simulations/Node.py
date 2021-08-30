from dataclasses import dataclass
from typing import Dict, List, Optional, Protocol

import numpy as np


class NodeLike(Protocol):
    def __hash__(self) -> int:
        ...

    def __str__(self) -> str:
        ...

    def construct_from_parent(self):
        ...

    def construct_from_source(self):
        ...


class SourcePopulation(object):
    rng = np.random.default_rng()

    def __init__(
        self, allele_frequencies: Dict[str, List[float]], mean_coi: float, label: str
    ) -> None:
        self.allele_frequencies = allele_frequencies
        self.mean_coi = mean_coi
        self.label = label

    def init_founder(self, node: NodeLike) -> NodeLike:
        node.infection_time = SourcePopulation.rng.uniform(1, 1000)
        node.observation_time = node.infection_time + node.infection_duration
        coi = SourcePopulation.rng.poisson(self.mean_coi) + 1
        for locus, dist in self.allele_frequencies.items():
            alleles = np.clip(SourcePopulation.rng.multinomial(coi, dist), 0, 1)
            node.true_alleles[locus] = alleles


class SimpleNode(object):
    _id_counter = 0
    rng = np.random.default_rng()

    def __init__(
        self,
        false_positive_rate: float,
        false_negative_rate: float,
        mutation_rate: float,
        loss_rate: float,
        infection_duration_shape: float,
        infection_duration_scale: float,
        parent: Optional["SimpleNode"] = None,
        source: Optional[SourcePopulation] = None,
        label: str = None,
    ) -> None:
        SimpleNode._id_counter += 1
        self.id = SimpleNode._id_counter
        self.label = label or str(self.id)
        self.false_positive_rate = false_positive_rate
        self.false_negative_rate = false_negative_rate
        self.mutation_rate = mutation_rate
        self.loss_rate = loss_rate
        self.infection_duration = SimpleNode.rng.gamma(
            infection_duration_shape, infection_duration_scale
        )
        self.parent = parent
        self.source = source

        self.infection_time = 0
        self.observation_time = 0

        self.true_alleles: Dict[str, List[float]] = {}
        self.observed_alleles: Dict[str, List[float]] = {}

        if parent:
            self.construct_from_parent()
        else:
            self.construct_from_source()

    def __hash__(self) -> int:
        return self.id

    def __str__(self) -> str:
        return self.label

    def __repr__(self) -> str:
        return f"<Node: {str(self)}>"

    def construct_from_parent(self):
        self.parent.init_child(self)
        self.generate_observed_alleles()

    def construct_from_source(self):
        self.source.init_founder(self)
        self.generate_observed_alleles()

    def generate_observed_alleles(self):
        for locus, alleles in self.true_alleles.items():
            observed_alleles = []
            for allele in alleles:
                if allele == 1:
                    observed_alleles.append(
                        SimpleNode.rng.binomial(1, 1 - self.false_negative_rate)
                    )
                else:
                    observed_alleles.append(
                        SimpleNode.rng.binomial(1, self.false_positive_rate)
                    )
            self.observed_alleles[locus] = np.array(observed_alleles)

    def init_child(self, child: "SimpleNode"):
        child.infection_time = SimpleNode.rng.uniform(
            self.infection_time, self.observation_time
        )
        child.observation_time = child.infection_time + child.infection_duration
        for locus, alleles in self.true_alleles.items():
            child_alleles = []
            total_alleles = sum(alleles)
            for allele in alleles:
                if allele == 1:
                    if total_alleles > 1:
                        child_alleles.append(
                            SimpleNode.rng.binomial(1, 1 - self.loss_rate)
                        )
                        total_alleles -= 1 - child_alleles[-1]
                    else:
                        child_alleles.append(allele)
                elif allele == 0:
                    child_alleles.append(SimpleNode.rng.binomial(1, self.mutation_rate))
            child.true_alleles[locus] = np.array(child_alleles)


def calculate_distance(node1, node2, loci=None) -> int:
    d = 0
    if not loci:
        loci = node1.observed_alleles.keys()
    for locus in loci:
        node1_a = node1.observed_alleles[locus]
        node2_a = node2.observed_alleles[locus]
        d += np.count_nonzero(node1_a != node2_a)
    return d
