import json
import numpy as np
from collections import defaultdict
from simulations.Node import SimpleNode, SourcePopulation, calculate_distance
from simulations.TransmissionNetwork import TransmissionNetwork, subsample_network
from typing import Dict


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


def parse_nodes(nodes, allowed_parents, loci=None):
    if not loci:
        loci = nodes[0].true_alleles.keys()
    nodes = [
        {
            "id": node.label,
            "latent_genotype": [
                {"locus": locus, "genotype": "".join(map(str, list(genotype)))}
                for (locus, genotype) in node.true_alleles.items()
                if locus in loci
            ],
            "observed_genotype": [
                {"locus": locus, "genotype": "".join(map(str, list(genotype)))}
                for (locus, genotype) in node.observed_alleles.items()
                if locus in loci
            ],
            "infection_duration": node.infection_duration,
            "infection_time": node.infection_time,
            "observation_time": node.observation_time,
            "allowed_parents": [p[0].label for p in allowed_parents[node]],
        }
        for node in nodes
    ]
    return nodes


def parse_network(edges):
    parents = edges.keys()
    network = []
    for parent in parents:
        for child in edges[parent]:
            network.append({"from": parent.label, "to": child.label})
    return network


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


def generate_allowed_parents(nodes, max_allowed_parents: int, loci=None) -> Dict:
    pairwise_distances = calculate_distance_matrix(nodes, loci)
    allowed_parents = {}
    for node in nodes:
        dists = pairwise_distances[node]
        dists.sort(key=lambda x: x[1])
        allowed_parents[node] = dists[0:max_allowed_parents]
    return allowed_parents


def calc_mean_bf(edge_sets):
    return np.mean([len(_) for _ in edge_sets])


rng = np.random.default_rng()

allele_frequencies = {}
for k in range(100):
    allele_frequencies[f"L{k + 1}"] = rng.dirichlet(np.array([1] * 10))

source_pop = SourcePopulation(
    allele_frequencies=allele_frequencies, mean_coi=4, label="S"
)

t1 = TransmissionNetwork(
    r0=0.825,
    node_constructor=SimpleNode,
    offspring_sampler=lambda x: rng.poisson(x),
    source_population=source_pop,
    num_founders=25,
)

# disallowed_parents = generate_allowed_parents(t1.nodes, 25)

# out = {
#     "loci": parse_loci(allele_frequencies),
#     "nodes": parse_nodes(t1, disallowed_parents),
#     "network": parse_network(t1),
# }

# with open(
#     "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/nodes11.json",
#     "w",
# ) as f:
#     json.dump(out, f, separators=(",", ":"))


# obs_50 = [subsample_network(t1, 0.5) for _ in range(10000)]
# obs_25 = [subsample_network(t1, 0.25) for _ in range(10000)]
# obs_125 = [subsample_network(t1, 0.125) for _ in range(10000)]
# obs_0625 = [subsample_network(t1, 0.0625) for _ in range(10000)]

# obs_90 = subsample_network(t1, 0.9)
# obs_75 = subsample_network(t1, 0.75)
# obs_50 = subsample_network(t1, 0.5)
# obs_25 = subsample_network(t1, 0.25)
# obs_125 = subsample_network(t1, 0.125)
# obs_0625 = subsample_network(t1, 0.0625)

# mean_bfs_50 = [calc_mean_bf(net[1].values()) for net in obs_50]
# mean_bfs_25 = [calc_mean_bf(net[1].values()) for net in obs_25]
# mean_bfs_125 = [calc_mean_bf(net[1].values()) for net in obs_125]
# mean_bfs_0625 = [calc_mean_bf(net[1].values()) for net in obs_0625]

# sub_network = subsample_network(t1, 0.5)


with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/superinf_2022_07_27/full_nodes.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(t1.nodes, generate_allowed_parents(t1.nodes, 25)),
            "network": parse_network(t1.edges),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/obs_exp/nodes_90.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(obs_90[0], generate_allowed_parents(obs_90[0], 25)),
            "network": parse_network(obs_90[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/obs_exp/nodes_75.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(obs_75[0], generate_allowed_parents(obs_75[0], 25)),
            "network": parse_network(obs_75[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/obs_exp/nodes_50.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(obs_50[0], generate_allowed_parents(obs_50[0], 25)),
            "network": parse_network(obs_50[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/obs_exp/nodes_25.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(obs_25[0], generate_allowed_parents(obs_25[0], 25)),
            "network": parse_network(obs_25[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/obs_exp/nodes_125.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies),
            "nodes": parse_nodes(
                obs_125[0], generate_allowed_parents(obs_125[0], 25)
            ),
            "network": parse_network(obs_125[1]),
        },
        f,
        separators=(",", ":"),
    )

# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
# ********************************************************************************

loci = list(allele_frequencies.keys())[0:25]

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/full_nodes.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(t1.nodes, generate_allowed_parents(t1.nodes, 25, loci), loci),
            "network": parse_network(t1.edges),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/nodes_90.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(obs_90[0], generate_allowed_parents(obs_90[0], 25, loci), loci),
            "network": parse_network(obs_90[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/nodes_75.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(obs_75[0], generate_allowed_parents(obs_75[0], 25, loci), loci),
            "network": parse_network(obs_75[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/nodes_50.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(obs_50[0], generate_allowed_parents(obs_50[0], 25, loci), loci),
            "network": parse_network(obs_50[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/nodes_25.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(obs_25[0], generate_allowed_parents(obs_25[0], 25, loci), loci),
            "network": parse_network(obs_25[1]),
        },
        f,
        separators=(",", ":"),
    )

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/locus_subset_exp/nodes_125.json",
        "w",
) as f:
    json.dump(
        {
            "loci": parse_loci(allele_frequencies, loci),
            "nodes": parse_nodes(
                obs_125[0], generate_allowed_parents(obs_125[0], 25, loci), loci
            ),
            "network": parse_network(obs_125[1]),
        },
        f,
        separators=(",", ":"),
    )
