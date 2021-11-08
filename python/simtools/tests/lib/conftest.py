import pytest

from simtools.lib.Network import Network, Node


@pytest.fixture(scope="module")
def dag():
    a = Node("a")
    b = Node("b")
    c = Node("c")
    d = Node("d")

    dag = Network()
    dag.add_edge(a.label, b.label)
    dag.add_edge(b.label, c.label)
    dag.add_edge(c.label, d.label)

    return dag
