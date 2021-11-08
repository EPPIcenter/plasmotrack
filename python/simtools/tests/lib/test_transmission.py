from simtools.lib.Transmission import (Infections,
                                       SimpleTransmissionProcess,
                                       SimpleSourceTransmissionProcess,
                                       generate_genetics_on_dag)


def test_transmission_infections(dag):
    stp = SimpleSourceTransmissionProcess(
        5, {"L1": [0.1, 0.2, 0.3, 0.4], "L2": [1 / 10] * 10, "L3": [1 / 5] * 5}, 25
    )

    gtp = SimpleTransmissionProcess(0.01, 25)

    generate_genetics_on_dag(dag, stp, gtp)
    assert len(Infections.get_infection("a").strains[0]["L1"]) == 4
    assert len(Infections.get_infection("b").strains[0]["L1"]) == 4
    assert len(Infections.get_infection("c").strains[0]["L1"]) == 4
    assert len(Infections.get_infection("d").strains[0]["L1"]) == 4
