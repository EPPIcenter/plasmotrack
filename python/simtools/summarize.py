from collections import defaultdict

import igraph
from matplotlib import pyplot as plt
import os
import os.path
import csv

output_dir = "/data/mmurphy/analysis/tnets/superinf_2022_09_22/output"

parent_set_files = os.listdir(os.path.join(output_dir, "parent_sets"))


edges = {}
processed_edges = []
for parent_set_file in parent_set_files:
    child = parent_set_file.split(".")[0].replace("_ps", "")
    edges[child] = defaultdict(float)
    count = 0
    with open(os.path.join(output_dir, "parent_sets", parent_set_file), "r") as f:
        r = csv.reader(f)
        next(r)
        while True:
            line = next(r)
            if len(line) > 1:
                if line[2] == "500":
                    break
        for line in r:
            if len(line) > 1:
                parent_set = line[0][1:-1].split(";")
                prob = float(line[1])
                for parent in parent_set:
                    edges[child][parent] += prob
            else:
                count += 1
    for parent in edges[child].keys():
        processed_edges.append(
            {"from": parent, "to": child, "weight": edges[child][parent] / count}
        )


graph = igraph.Graph.DictList(
    vertices=[{"id": node} for node in edges.keys() if "copy" not in node],
    edges=[edge for edge in processed_edges if "copy" not in edge["from"]],
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
    edge_color=[f'rgba(1,1,1,{edge["weight"]})' for edge in graph.es],
    # edge_weight=graph.es["weight"],
    vertex_label=graph.vs["id"],
)
plt.show()
