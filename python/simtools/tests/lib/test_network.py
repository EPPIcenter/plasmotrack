from simtools.lib.Network import Network, Node


def test_node_ids_unique():
    a = Node("a")
    b = Node("b")
    assert a.label != b.label


def test_nodes_tracked():
    a = Node("a")
    b = Node("b")

    assert a == Node.get_node("a")
    assert b == Node.get_node("b")
    assert a.label != b.label


def test_network_created():
    a = Node("a")
    b = Node("b")
    c = Node("c")
    d = Node("d")

    dag = Network()
    dag.add_edge(a.label, b.label)
    dag.add_edge(b.label, c.label)
    dag.add_edge(c.label, d.label)

    assert not dag.parent_sets[a.label]
    assert len(dag.parent_sets[b.label]) == 1
    assert len(dag.parent_sets[c.label]) == 1
    assert len(dag.parent_sets[d.label]) == 1
    assert len(dag.get_roots()) == 1
