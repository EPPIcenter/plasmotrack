import numpy as np
from copy import deepcopy
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
    def __init__(self, loss_rate: float, seed: int = 0) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)
        self.loss_rate = loss_rate

    def apply(self, infections: List[Infection]) -> Infection:
        latent_infection: Infection = [
            deepcopy(strain) for infection in infections for strain in infection
        ]

        # for strain in all_strains:
        #     if self.rng.binomial(1, 1 - self.loss_rate):
        #         latent_infection.append(strain)

        gametocytes = random.choices(latent_infection, k=np.random.poisson(5) + 1)
        num_oocyst = np.random.poisson(5) + 1
        strains = []
        for _ in range(num_oocyst):
            strains.append(recombine_strains(gametocytes))

        num_spor = np.random.poisson(3) + 1
        transmitted_infection = random.choices(strains, k=num_spor)

        return transmitted_infection


def recombine_strains(infection: Infection):
    s1 = random.choice(infection)
    s2 = random.choice(infection)
    s_out = deepcopy(s1)
    for locus in s_out.keys():
        if random.random() > 0.5:
            s_out[locus] = deepcopy(s2[locus])

    return s_out


def generate_genetics_on_dag(
    dag: Network, stp: SourceTransmissionProcess, tp: TransmissionProcess
) -> None:
    to_process = dag.get_topological_sorting()

    for node in to_process:
        parent_set = dag.get_parents(node)
        if not parent_set:
            gp = stp.apply()
            Infections(node, gp)
        else:
            parent_set_infections = [
                Infections.get_infection(p).strains for p in parent_set
            ]
            inf = tp.apply(parent_set_infections)
            Infections(node, inf)
