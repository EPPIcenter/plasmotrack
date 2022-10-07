import numpy as np
from collections import defaultdict
from simtools.lib.Detection import InfectionDuration, InfectionTime
from simtools.lib.Observation import ObservedGenetics
from simtools.lib.Transmission import Infections


def calculate_distance(node1_label, node2_label, loci=None) -> int:
    node1_genetics = ObservedGenetics.get_observed_genetics(node1_label)
    node2_genetics = ObservedGenetics.get_observed_genetics(node2_label)
    d = 0
    if not loci:
        loci = node1_genetics.genotypes.keys()
    for locus in loci:
        node1_a = node1_genetics.genotypes[locus]
        node2_a = node2_genetics.genotypes[locus]
        d += np.count_nonzero(np.logical_and(node1_a, node2_a))
    return d


def parse_network(dag):
    parents = dag.parent_sets.keys()
    network = []
    for parent in parents:
        for child in dag.get_children(parent):
            network.append({"from": parent, "to": child})
    return network


def parse_loci(allele_frequencies, allowed_loci=None):
    if not allowed_loci:
        allowed_loci = allele_frequencies.keys()
    loci = [
        {
            "locus": locus,
            "allele_freqs": list(allele_freqs),
            "num_alleles": len(allele_freqs),
        }
        for (locus, allele_freqs) in allele_frequencies.items()
        if locus in allowed_loci
    ]
    return loci


def parse_nodes(node_labels, allowed_parents, loci=None):
    if not loci:
        loci = ObservedGenetics.get_observed_genetics(node_labels[0]).genotypes.keys()
    nodes = [
        {
            "id": node,
            "latent_genotype": [
                [
                    {"locus": locus, "genotype": genotype.tolist()}
                    for locus, genotype in strain.items()
                    if locus in loci
                ]
                for strain in Infections.get_infection(node).strains
            ],
            "flat_latent_genotype": [
                {
                    "locus": locus,
                    "genotype": flatten_genotypes(
                        Infections.get_infection(node).strains, locus
                    ),
                }
                for locus in loci
            ],
            "observed_genotype": [
                {"locus": locus, "genotype": "".join(map(str, list(genotype)))}
                for locus, genotype in ObservedGenetics.get_observed_genetics(
                    node
                ).genotypes.items()
                if locus in loci
            ],
            "infection_duration": InfectionDuration.get_infection_duration(
                node
            ).duration,
            "infection_time": InfectionTime.get_infection_time(node).time,
            "observation_time": InfectionTime.get_infection_time(node).time
            + InfectionDuration.get_infection_duration(node).duration,
            "allowed_parents": [p[0] for p in allowed_parents[node]],
        }
        for node in node_labels
    ]
    return nodes


def flatten_genotypes(strains, locus):
    genotypes = [strain[locus] for strain in strains]
    flat_genotype = np.sum(genotypes, axis=0)
    flat_genotype = flat_genotype > 0
    flat_genotype = "".join(map(str, list(flat_genotype.astype(int))))
    return flat_genotype


def calculate_distance_matrix(nodes, loci=None):
    pairwise_distances = defaultdict(list)
    total_nodes = len(nodes)

    for i in range(total_nodes):
        for j in range(i + 1, total_nodes):
            node_i = nodes[i]
            node_j = nodes[j]
            d = calculate_distance(node_i, node_j, loci)
            pairwise_distances[node_i].append((node_j, d))
            pairwise_distances[node_j].append((node_i, d))
    return pairwise_distances


def generate_allowed_parents(nodes, max_allowed_parents: int, network, loci=None):
    pairwise_distances = calculate_distance_matrix(nodes, loci)
    allowed_parents = {}
    for node in nodes:
        dists = pairwise_distances[node]
        dists.sort(key=lambda x: -x[1])
        subset = dists[:max_allowed_parents]
        true_parents = [p["from"] for p in network if p["to"] == node]
        true_parents = [_ for _ in dists if _[0] in true_parents]
        print(true_parents)
        allowed_parents[node] = list(set(subset + true_parents))

        print(node, allowed_parents[node])
    return allowed_parents
