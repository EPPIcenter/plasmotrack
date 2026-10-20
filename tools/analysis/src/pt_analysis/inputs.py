"""The JSON files an Analysis names besides its stores: the engine input and the simulation's truth, plain or
gzipped (`.gz`).

Engine input: `loci` [{locus, num_alleles, allele_freqs}] and `nodes` [{id, observation_time, symptomatic,
observed_genotype: [{locus, genotype}], allowed_parents}]. A genotype string's character k is Allele k, `1` present;
an all-zero, empty or missing genotype is unobserved. Truth: `network` [{src, dest}] and `nodes` with
`flat_latent_genotype` [{locus, genotype}], `infection_time`, `infection_duration` and `observation_time`.
"""

import gzip
import json
from functools import cache
from pathlib import Path


@cache
def _read(path, mtime):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return json.load(f)


def read_json(path):
    """The JSON file at `path`, gunzipped when its name ends in `.gz`; cached until the file changes."""
    return _read(str(path), Path(path).stat().st_mtime_ns)


def nodes(path):
    """{Infection id (as a string): node} of an input or truth file."""
    return {str(node["id"]): node for node in read_json(path)["nodes"]}


def candidate_parents(path):
    """{Infection id: [its `allowed_parents` that are Infections of the input, as strings]} of an input file."""
    found = nodes(path)
    return {child: [str(p) for p in node.get("allowed_parents") or [] if str(p) in found]
            for child, node in found.items()}


def loci(path):
    """The Locus names of an input file, in its order."""
    return [locus["locus"] for locus in read_json(path)["loci"]]


def genotypes(entries):
    """{locus: genotype string} from a node's `observed_genotype` or `flat_latent_genotype` list, leaving out
    unobserved Loci (empty or all-zero strings)."""
    return {e["locus"]: e["genotype"] for e in entries or [] if "1" in e.get("genotype", "")}
