from simtools.lib.Detection import (InfectionDuration, InfectionTime,
                                    PassiveDetectionProcess,
                                    UniformInfectionTimeProcess,
                                    UniformSourceInfectionTimeProcess,
                                    generate_infection_times_on_dag)


def test_transmission_infection_time(dag):
    sitp = UniformSourceInfectionTimeProcess(-1, 1000, 25)
    itp = UniformInfectionTimeProcess(9, 25)
    dtp = PassiveDetectionProcess(9, 10, 25)

    generate_infection_times_on_dag(dag, sitp, itp, dtp)
    assert InfectionDuration.get_infection_duration("a").duration > 0
    assert InfectionDuration.get_infection_duration("b").duration > 0
    assert InfectionDuration.get_infection_duration("c").duration > 0
    assert InfectionDuration.get_infection_duration("d").duration > 0
    assert (
        InfectionTime.get_infection_time("a").time
        < InfectionTime.get_infection_time("b").time
    )
    assert (
        InfectionTime.get_infection_time("b").time
        < InfectionTime.get_infection_time("c").time
    )
    assert (
        InfectionTime.get_infection_time("c").time
        < InfectionTime.get_infection_time("d").time
    )
