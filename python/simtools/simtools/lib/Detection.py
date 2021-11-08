from typing import Dict, Protocol

import numpy as np

from simtools.lib.Network import Network


class InfectionTime:
    _nodes: Dict[str, "InfectionTime"] = {}

    def __init__(self, label: str, time: float) -> None:
        self.label = label
        self.time = time
        InfectionTime._nodes[self.label] = self

    @classmethod
    def get_infection_time(cls, label: str) -> "InfectionTime":
        return cls._nodes[label]


class InfectionDuration:
    _nodes: Dict[str, "InfectionDuration"] = {}

    def __init__(self, label: str, duration: float) -> None:
        self.label = label
        self.duration = duration
        InfectionDuration._nodes[self.label] = self

    @classmethod
    def get_infection_duration(cls, label: str) -> "InfectionDuration":
        return cls._nodes[label]


class DetectionProcess(Protocol):
    def apply(self, infection_epidemiology: Dict[str, any] = None) -> float:
        pass


class PassiveDetectionProcess(DetectionProcess):
    def __init__(self, shape: float, scale: float, seed: int = 0) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)
        self.shape = shape
        self.scale = scale

    def apply(self, _: Dict[str, any] = None) -> float:
        return self.rng.gamma(self.shape, self.scale)


class InfectionTimeProcess(Protocol):
    """
    Given the parent time of infection generates the time of infection for a child infection
    """

    def apply(self, parent_time_of_infection: float, parent_time_of_detection: float):
        pass


class SourceInfectionTimeProcess(Protocol):
    def apply(self) -> float:
        pass


class UniformInfectionTimeProcess(InfectionTimeProcess):
    """
    Uniform probability of transmission over duration of infection
    """

    def __init__(self, min_duration: float, seed: int = -1) -> None:
        self.seed = seed
        self.min_duration = min_duration
        self.rng = np.random.default_rng(seed=seed)

    def apply(
        self, parent_time_of_infection: float, parent_infection_duration: float
    ) -> float:
        return parent_time_of_infection + self.rng.uniform(
            self.min_duration, parent_infection_duration
        )


class UniformSourceInfectionTimeProcess(SourceInfectionTimeProcess):
    def __init__(self, min_time: float, max_time: float, seed: int = -1) -> None:
        self.seed = seed
        self.min_time = min_time
        self.max_time = max_time
        self.rng = np.random.default_rng(seed=seed)

    def apply(self) -> float:
        return self.rng.uniform(self.min_time, self.max_time)


def generate_infection_times_on_dag(
    dag: Network,
    sitp: SourceInfectionTimeProcess,
    itp: InfectionTimeProcess,
    dtp: DetectionProcess,
) -> None:
    source_nodes = dag.get_roots()
    to_process = dag.get_topological_sorting()
    for node in source_nodes:
        to_process += dag.child_sets[node]
    for node in to_process:
        parent_set = dag.parent_sets[node]
        if not parent_set: 
            InfectionTime(node, sitp.apply())
            InfectionDuration(node, dtp.apply())
        else:
            parent_set_infection_times = [
                InfectionTime.get_infection_time(p) for p in parent_set
            ]
            parent_set_infection_durations = [
                InfectionDuration.get_infection_duration(p) for p in parent_set
            ]
            infection_time = max(
                [
                    itp.apply(i.time, d.duration)
                    for i, d in zip(
                        parent_set_infection_times, parent_set_infection_durations
                    )
                ]
            )
            InfectionTime(node, infection_time)
            InfectionDuration(node, dtp.apply())
