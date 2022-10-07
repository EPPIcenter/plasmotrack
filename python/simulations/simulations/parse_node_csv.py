import csv
import json
import numpy as np
from collections import defaultdict


def parse_node_csv(file_path):
    nodes = []
    ignore = ["Type", "History", "Lat", "Lon", "gpslat", "gpslon"]
    with open(file_path, "r") as f:
        r = csv.DictReader(f)
        i = 1
        for node in r:
            node_entry = {"id": str(i), "observation_time": int(node.pop("Time"))}
            # node_entry = {"id": str(i), "observation_time": int(node.pop("Date"))}
            node_entry["observed_genotype"] = [
                {"locus": locus, "genotype": node[locus]} for locus in node.keys() if locus not in ignore
            ]
            nodes.append(node_entry)
            i += 1
    return nodes


def calculate_distance(node1, node2):
    d = 0
    loci1 = set([el["locus"] for el in node1["observed_genotype"]])
    loci2 = set([el["locus"] for el in node2["observed_genotype"]])
    loci = loci1.intersection(loci2)
    for locus in loci:
        node1_a = [
            loc["genotype"]
            for loc in node1["observed_genotype"]
            if loc["locus"] == locus
        ][0]
        node2_a = [
            loc["genotype"]
            for loc in node2["observed_genotype"]
            if loc["locus"] == locus
        ][0]
        d += np.count_nonzero(node1_a != node2_a)
    return d


def calculate_distance_matrix(nodes):
    pairwise_distances = defaultdict(list)
    total_nodes = len(nodes)

    for i in range(total_nodes):
        print(i)
        for j in range(i + 1, total_nodes):
            node_i = nodes[i]
            node_j = nodes[j]
            d = calculate_distance(node_i, node_j)
            pairwise_distances[node_i["id"]].append((node_j["id"], d))
            pairwise_distances[node_j["id"]].append((node_i["id"], d))
    return pairwise_distances


def generate_allowed_parents(nodes, max_allowed_parents):
    pairwise_distances = calculate_distance_matrix(nodes)
    allowed_parents = {}
    for node in nodes:
        dists = pairwise_distances[node["id"]]
        dists.sort(key=lambda x: x[1])
        allowed_parents[node["id"]] = dists[0:max_allowed_parents]
    return allowed_parents


parsed_nodes = parse_node_csv(
    "/home/mmurphy/analysis/Zanzibar2021/Data/Microsat/Znz_nodes_data20.csv"
)
allowed_parents = generate_allowed_parents(parsed_nodes, 15)

for node in parsed_nodes:
    node["allowed_parents"] = [str(x[0]) for x in allowed_parents[node["id"]]]

loci = [
    {"locus": l["locus"], "num_alleles": len(l["genotype"])}
    for l in parsed_nodes[0]["observed_genotype"]
]

with open(
        "/home/mmurphy/Dropbox/transmission_nets_test_files/resources/JSON/zanzibar/nodes.json",
        "w",
) as f:
    json.dump({"loci": loci, "nodes": parsed_nodes}, f, separators=(",", ":"))
