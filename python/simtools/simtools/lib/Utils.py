import numpy as np

from typing import Dict
from simtools.lib.Observation import Genotype


def calculate_distance(a: Dict[str, Genotype], b: Dict[str, Genotype]):
    d = 0
    common_loci = set.intersection(set(a.keys()), set(b.keys()))
    for locus in common_loci:
        d += np.count_nonzero(np.array(a[locus]) != np.array(b[locus]))
    return d
