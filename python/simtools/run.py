from simtools.lib.Network import Node, Network
from simtools.lib.Detection import (InfectionDuration, InfectionTime,
                                    PassiveDetectionProcess,
                                    UniformInfectionTimeProcess,
                                    UniformSourceInfectionTimeProcess,
                                    generate_infection_times_on_dag)
from simtools.lib.Transmission import SimpleTransmissionProcess, SimpleSourceTransmissionProcess, Infections, generate_genetics_on_dag
from simtools.lib.Observation import SimpleGeneticsObservationProcess, generate_observed_genetics, ObservedGenetics

import random
import numpy as np
from collections import defaultdict
import json


def calculate_distance(node1_label, node2_label, loci=None) -> int:
    node1_genetics = ObservedGenetics.get_observed_genetics(node1_label)
    node2_genetics = ObservedGenetics.get_observed_genetics(node2_label)
    d = 0
    if not loci:
        loci = node1_genetics.genotypes.keys()
    for locus in loci:
        node1_a = node1_genetics.genotypes[locus]
        node2_a = node2_genetics.genotypes[locus]
        d += np.count_nonzero(node1_a != node2_a)
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
                [{"locus": locus, "genotype": genotype.tolist()} for locus, genotype in strain.items() if locus in loci] for strain in Infections.get_infection(node).strains
            ],
            "observed_genotype": [
                {"locus": locus, "genotype" : "".join(map(str, list(genotype)))} for locus, genotype in ObservedGenetics.get_observed_genetics(node).genotypes.items() if locus in loci 
            ],
            "infection_duration": InfectionDuration.get_infection_duration(node).duration,
            "infection_time": InfectionTime.get_infection_time(node).time,
            "observation_time": InfectionTime.get_infection_time(node).time + InfectionDuration.get_infection_duration(node).duration,
            "allowed_parents": [p[0] for p in allowed_parents[node]],
        }
        for node in node_labels
    ]
    return nodes

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


def generate_allowed_parents(nodes, max_allowed_parents: int, loci=None):
    pairwise_distances = calculate_distance_matrix(nodes, loci)
    allowed_parents = {}
    for node in nodes:
        dists = pairwise_distances[node]
        dists.sort(key=lambda x: x[1])
        allowed_parents[node] = dists[0:max_allowed_parents]
    return allowed_parents



dag = Network()

idx = 1

outbreaks = []

for i in range(idx,11):
    node_id = str(i)
    node = Node(node_id)
    idx += 1
    outbreaks.append([node.label])

for i in range(idx, 100):
    node_id = str(i)
    outbreak_choice = random.choice(range(len(outbreaks)))
    # num_parents = 1 + 1 * (random.random() > .95)
    # num_parents = random.choice([1,2])
    num_parents = 1
    parents = list(set(random.choices(outbreaks[outbreak_choice], k=num_parents)))

    node = Node(node_id)
    outbreaks[outbreak_choice].append(node.label)
    for parent_id in parents:
        p_node = Node.get_node(parent_id)
        dag.add_edge(p_node.label, node.label)
    
loci = {"L" + str(i) : [1/15] * 15 for i in range(1, 101)}


sitp = UniformSourceInfectionTimeProcess(0, 1000, 25)
itp = UniformInfectionTimeProcess(9, 25)
dtp = PassiveDetectionProcess(9, 10, 25)
gtp = SimpleTransmissionProcess(0.1, 25)
stp = SimpleSourceTransmissionProcess(5, loci)
op = SimpleGeneticsObservationProcess(.01, .01)

generate_genetics_on_dag(dag, stp, gtp)
generate_infection_times_on_dag(dag, sitp, itp, dtp)
generate_observed_genetics(dag, op)


ap = generate_allowed_parents(list(Node._nodes.keys()), 15)
nodes = list(Node._nodes.keys())


with open(
    "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/superinf/full_nodes.json",
    "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(loci),
            "nodes": parse_nodes(nodes, ap),
            "network": parse_network(dag),
        },
        f,
        separators=(",", ":"),
    )