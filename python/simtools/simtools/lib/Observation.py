"""
The observation module implements the observation process that samples undergo. The observed state of a sample is a noisy
representation of some latent state including genetics, and potentially the time of observation and geography.
"""
from typing import Dict, List, Protocol

import numpy as np
from simtools.lib.Network import Network, Node

from simtools.lib.Transmission import Infection, Infections, Locus

Genotype = Dict[str, Locus]


class ObservedGenetics:
    _nodes: Dict[str, "ObservedGenetics"] = {}

    def __init__(self, label: str, genotypes: Dict[str, Genotype]) -> None:
        self.label = label
        self.genotypes = genotypes
        ObservedGenetics._nodes[self.label] = self

    @classmethod
    def get_observed_genetics(cls, label: str) -> "ObservedGenetics":
        return cls._nodes[label]


class GeneticsObservationProcess(Protocol):
    def apply(self, infection: Infection) -> Genotype:
        pass


class SimpleGeneticsObservationProcess(GeneticsObservationProcess):
    def __init__(self, fpr: float, fnr: float, seed: int = 0) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)
        self.fpr = fpr
        self.fnr = fnr

    def _flatten_locus(self, strains: List[Locus]) -> Locus:
        out = [0] * len(strains[0]) 
        for genotype in strains:
            for idx, allele in enumerate(genotype.tolist()):
                out[idx] = out[idx] or allele
        return out

    def _flatten_strains(self, infection: Infection) -> Genotype:
        loci = infection[0].keys()
        out = {
            locus: self._flatten_locus([strain[locus] for strain in infection])
            for locus in loci
        }
        return out

    def apply(self, infection: Infection) -> Genotype:
        flat_genotype = self._flatten_strains(infection)
        out: Genotype = {}
        for label, locus in flat_genotype.items():
            obs_locus = []
            for allele in locus:
                if allele == 1:
                    obs_locus.append(self.rng.binomial(1, 1 - self.fnr))
                else:
                    obs_locus.append(self.rng.binomial(1, self.fpr))
            out[label] = obs_locus
        return out


def generate_observed_genetics(dag, op: GeneticsObservationProcess):
    to_process = dag.get_topological_sorting()
    for label in to_process:
        infection = Infections.get_infection(label)
        obs_genetics = op.apply(infection.strains)
        ObservedGenetics(label, obs_genetics)