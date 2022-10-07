from simtools.lib.Observation import ObservedGenetics
from simtools.lib.Utils import calculate_distance


def test_distance_metric():
    a = ObservedGenetics('a', {
        'L1': [1, 0, 0, 1],
        'L2': [0, 0, 0, 1]
    })
    b = ObservedGenetics('b', {
        'L1': [1, 1, 0, 1],
        'L2': [1, 0, 0, 1]
    })
    c = ObservedGenetics('c', {
        'L1': [1, 1, 1, 1],
        'L3': [1, 0, 0, 1]
    })

    assert calculate_distance(a.genotypes, b.genotypes) == 2
    assert calculate_distance(a.genotypes, c.genotypes) == 2
    assert calculate_distance(b.genotypes, c.genotypes) == 1
