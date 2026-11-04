"""Times of Infections: observation times from the engine input, the lag between two Infections, and the posterior
infection time pooled over an Analysis's Chains."""

from functools import reduce

from pt_analysis import inputs
from pt_analysis.derive import table
from pt_analysis.sql import Scope


def observation_times(analysis):
    """{Infection id: observation time} from the Analysis's input, None for a node without one."""
    return {i: node.get("observation_time") for i, node in inputs.nodes(analysis.input).items()}


def lag(times, parent, child):
    """The child's observation time minus the parent's, from `observation_times`; None when either has none."""
    if times.get(parent) is None or times.get(child) is None:
        return None
    return times[child] - times[parent]


def infection_times(analysis, infections=None):
    """{Infection: {lower, median, upper}}: the posterior infection time, the observation time minus the 97.5%, 50%
    and 2.5% quantiles (type 7, `quantile_cont`) of `infection_duration` pooled over every Draw of every Chain. For
    `infections` only when given; an Infection without Draws or without an observation time is left out."""
    times = observation_times(analysis)
    durations = reduce(lambda a, b: a.union(b), (table(store, "infections").project("infection, infection_duration")
                                                 for store in analysis.runs.values()))
    with Scope(next(iter(analysis.runs.values())).connection) as s:
        rows = s.rows(f"""
            SELECT infection, quantile_cont(infection_duration, 0.975), quantile_cont(infection_duration, 0.5),
                   quantile_cont(infection_duration, 0.025)
            FROM {s.name(durations)} WHERE $ids IS NULL OR list_contains($ids, infection)
            GROUP BY infection""", {"ids": None if infections is None else list(infections)})
    return {infection: dict(zip(("lower", "median", "upper"), (times[infection] - d for d in quantiles)))
            for infection, *quantiles in rows if times.get(infection) is not None}
