import json
import random
import os
import os.path

from simtools.lib.Detection import (
    PassiveDetectionProcess,
    UniformInfectionTimeProcess,
    UniformSourceInfectionTimeProcess,
    generate_infection_times_on_dag,
)
from simtools.lib.Network import Node, Network
from simtools.lib.Observation import (
    SimpleGeneticsObservationProcess1,
    generate_observed_genetics,
)
from simtools.lib.Transmission import (
    SimpleTransmissionProcess,
    SimpleSourceTransmissionProcess,
    generate_genetics_on_dag,
)
from simtools.lib.Utils import *

seed = 0

dag = Network()

num_outbreaks = 50
total_nodes = 200
num_alleles = 10
num_loci = 100
mean_coi = 3
loss_rate = 0.1
fpr = 0.01
fnr = 0.01


outbreaks = []
node_idx = 1
for i in range(0, num_outbreaks):
    outbreaks.append([])
    src_nodes_per_outbreak = random.choice([1, 2, 3])
    print(f"Outbreak {i} has {src_nodes_per_outbreak} source nodes")
    for _ in range(src_nodes_per_outbreak):
        node = Node(str(node_idx))
        node_idx += 1
        outbreaks[i].append(node.label)

print(outbreaks)

superinfection_prob = 0.2

for i in range(node_idx, total_nodes):
    # outbreak_choice = random.choice(range(num_outbreaks - 2))
    outbreak_choice = random.choice(range(25))
    num_parents = min(
        len(outbreaks[outbreak_choice]),
        1 + 1 * (random.random() > (1 - superinfection_prob)),
    )
    parents = list(set(random.choices(outbreaks[outbreak_choice], k=num_parents)))

    node = Node(str(i))
    outbreaks[outbreak_choice].append(node.label)
    for parent_id in parents:
        p_node = Node.get_node(parent_id)
        dag.add_edge(p_node.label, node.label)

print(outbreaks)


loci = {"L" + str(i): [1 / num_alleles] * num_alleles for i in range(1, num_loci + 1)}

sitp = UniformSourceInfectionTimeProcess(0, 1000, seed)
itp = UniformInfectionTimeProcess(9, seed)
dtp = PassiveDetectionProcess(9, 10, seed)
gtp = SimpleTransmissionProcess(loss_rate, seed)
stp = SimpleSourceTransmissionProcess(mean_coi, loci, seed)
op = SimpleGeneticsObservationProcess1(fpr, fnr, seed)

generate_genetics_on_dag(dag, stp, gtp)
generate_infection_times_on_dag(dag, sitp, itp, dtp)
generate_observed_genetics(dag, op)

network = parse_network(dag)
ap = generate_allowed_parents(list(Node._nodes.keys()), 10, network)
nodes = list(Node._nodes.keys())

output = "/data/mmurphy/analysis/tnets/superinf_2022_10_05/full_nodes.json"
os.makedirs(os.path.dirname(output), exist_ok=True)

with open(output, "w") as f:
    json.dump(
        {
            "loci": parse_loci(loci),
            "nodes": parse_nodes(nodes, ap),
            "network": parse_network(dag),
        },
        f,
        separators=(",", ":"),
    )
