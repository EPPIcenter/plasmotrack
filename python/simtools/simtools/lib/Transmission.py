from copy import deepcopy

import numpy as np

from simtools.lib.Network import Network

"""
Given the network topology, overlay elements of transmission. This includes genetics, time between infections, and infection geography
"""

from typing import Dict, List, Protocol
import random

Locus = List[int]  # Binary sequence indicating presence/absence of alleles
Strain = Dict[str, Locus]  # Collection of loci make up a strain
Infection = List[Strain]  # Collection of strains make up an infection

DiscreteDist = List[float]
AlleleFrequencies = Dict[str, DiscreteDist]


class Infections:
    _nodes: Dict[str, "Infections"] = {}

    def __init__(self, label: str, strains: Infection) -> None:
        self.label = label
        self.strains = strains
        Infections._nodes[self.label] = self

    @classmethod
    def get_infection(cls, label: str) -> "Infections":
        return cls._nodes[label]


class TransmissionProcess(Protocol):
    """
    Takes as input a list of infections (i.e. the parent set) and returns a new infection
    """

    def apply(self, infections: List[Infection]) -> Infection:
        pass


class SourceTransmissionProcess(Protocol):
    """
    Generates a new infection sampled from a source population.
    """

    def apply(self) -> Infection:
        pass


class SimpleSourceTransmissionProcess(SourceTransmissionProcess):
    def __init__(
        self, lam: float, allele_frequencies: AlleleFrequencies, seed: int = 0
    ) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)
        self.lam = lam
        self.allele_frequencies = allele_frequencies

    def apply(self) -> Infection:
        coi = self.rng.poisson(self.lam) + 1
        strains: Infection = []
        for _ in range(coi):
            strain: Strain = {
                locus: self.rng.multinomial(1, self.allele_frequencies[locus])
                for locus in self.allele_frequencies
            }
            strains.append(strain)
        return strains


class SimpleTransmissionProcess(TransmissionProcess):
    """
    No super infection
    """

    def __init__(self, loss_rate: float, seed: int = 0) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)
        self.loss_rate = loss_rate

    def apply(self, infections: List[Infection]) -> Infection:
        transmitted_infection: Infection = []
        all_strains = [deepcopy(strain) for infection in infections for strain in infection]
        transmitted_infection.append(all_strains.pop(random.choice(range(len(all_strains)))))
        for strain in all_strains:
            if self.rng.binomial(1, 1 - self.loss_rate):
                transmitted_infection.append(strain)
        # for infection in infections:
        #     for strain in infection:
        #         if self.rng.binomial(1, 1 - self.loss_rate):
        #             transmitted_infection.append(deepcopy(strain))
        return transmitted_infection


def generate_genetics_on_dag(
    dag: Network, stp: SourceTransmissionProcess, tp: TransmissionProcess
) -> None:
    to_process = dag.get_topological_sorting()
        
    for node in to_process:
        # to_process += dag.child_sets[node]
        parent_set = dag.parent_sets[node]
        if not parent_set:
            gp = stp.apply()
            Infections(node, gp)
        else:
            parent_set_infections = [
                Infections.get_infection(p).strains for p in parent_set
            ]
            inf = tp.apply(parent_set_infections)
            Infections(node, inf)
