import json
import numpy as np
import igraph
from matplotlib import pyplot as plt


def to_binary(genotype):
    return np.array([1 if allele == "1" else 0 for allele in genotype])


def calulate_similarity(node_a, node_b):
    similarity = 0
    pairs = zip(node_a, node_b)
    for pair in pairs:
        similarity += np.count_nonzero(
            np.logical_and(
                to_binary(pair[0]["genotype"]), to_binary(pair[1]["genotype"])
            )
        )
    return similarity


def calculate_difference(node_a, node_b):
    difference = 0
    pairs = zip(node_a, node_b)
    for pair in pairs:
        difference += np.count_nonzero(
            np.logical_xor(
                to_binary(pair[0]["genotype"]), to_binary(pair[1]["genotype"])
            )
        )
    return difference


input_file = "/data/mmurphy/analysis/tnets/superinf_2022_09_22/full_nodes.json"

with open(input_file, "r") as f:
    data = json.load(f)


nodes = data["nodes"]
network = data["network"]

node_dict = {node["id"]: node for node in nodes}

vertices = [
    {
        "id": node["id"],
        "observation_time": round(node["observation_time"]),
        "infection_time": round(node["infection_time"]),
        "infection_duration": round(node["infection_duration"]),
    }
    for node in nodes
]

for edge in network:
    from_node = node_dict[edge["from"]]
    to_node = node_dict[edge["to"]]
    true_similarity = calulate_similarity(
        from_node["flat_latent_genotype"], to_node["flat_latent_genotype"]
    )
    obs_similarity = calulate_similarity(
        from_node["observed_genotype"], to_node["observed_genotype"]
    )
    true_difference = calculate_difference(
        from_node["flat_latent_genotype"], to_node["flat_latent_genotype"]
    )
    obs_difference = calculate_difference(
        from_node["observed_genotype"], to_node["observed_genotype"]
    )
    edge["true_similarity"] = true_similarity
    edge["obs_similarity"] = obs_similarity
    edge["true_difference"] = true_difference
    edge["obs_difference"] = obs_difference

graph = igraph.Graph.DictList(
    vertices=vertices,
    edges=network,
    directed=True,
    vertex_name_attr="id",
    edge_foreign_keys=("from", "to"),
)

layout = graph.layout("fruchterman_reingold")
fig, ax = plt.subplots()
igraph.plot(
    graph,
    layout=layout,
    target=ax,
    # edge_label=graph.es["obs_similarity"],
    # edge_label=graph.es["obs_difference"],
    # edge_label=graph.es["true_difference"],
    # vertex_label=graph.vs["infection_time"],
    vertex_label=graph.vs["id"],
    vertex_size=0,
)
plt.show()
