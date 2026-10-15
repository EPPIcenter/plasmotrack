"""The sampler panel's data (U4), for one Run: are the samplers doing their jobs?

Blocks (`acceptance`, `swaps`) hold counts per Window, not per Draw, so they are read over both Phases, and each row
carries its Window's Phase and Sampling-iteration range. Run events come from `events`, whose payloads are JSON
objects (`doc/store-layout.md`).
"""

import json

from pt_analysis.store import BOTH_PHASES

# Run events marked on every time axis of the panel.
MARKED = ("non_finite_value", "interrupted", "burn_in_ended", "run_resumed", "run_ended")


def _windows(store):
    return {w.number: w for w in store.closed_windows()}


def _rate(accepted, total):
    return accepted / total if total else None


def acceptance(store):
    """One row per sampler per Window, over both Phases: window, phase, first_iteration, last_iteration, sampler,
    type, target, accepted, proposed and rate (accepted / proposed; None when nothing was proposed). Type and target
    come from `structure/samplers`. Ordered by type, target, then Window."""
    relation = store.table("acceptance", BOTH_PHASES)
    if relation is None:
        return []
    split = {s: (t, g) for s, t, g in store.structure("samplers").fetchall()}
    windows = _windows(store)
    rows = []
    for window, phase, name, accepted, proposed in relation.select(
            '"window", phase, sampler, accepted, proposed').fetchall():
        kind, target = split.get(name) or _split(name)
        w = windows[window]
        rows.append({"window": window, "phase": phase, "first_iteration": w.first_iteration,
                     "last_iteration": w.last_iteration, "sampler": name, "type": kind, "target": target,
                     "accepted": accepted, "proposed": proposed, "rate": _rate(accepted, proposed)})
    rows.sort(key=lambda r: (r["type"], r["target"], r["window"]))
    return rows


def _split(name):
    """A sampler id split at its first `/`, as `structure/samplers` does, for an id missing from that table."""
    kind, _, target = name.partition("/")
    return kind, target


def swaps(store):
    """One row per neighbouring pair of temperature slots per Window, over both Phases: window, phase,
    first_iteration, last_iteration, temperature_pair (k: slots k and k + 1), accepted, attempted and rate (None when
    nothing was attempted). Ordered by Window, then pair. [] for a Run with one chain, which has no `swaps`."""
    relation = store.table("swaps", BOTH_PHASES)
    if relation is None:
        return []
    windows = _windows(store)
    rows = []
    for window, phase, pair, accepted, attempted in relation.select(
            '"window", phase, temperature_pair, accepted, attempted').fetchall():
        w = windows[window]
        rows.append({"window": window, "phase": phase, "first_iteration": w.first_iteration,
                     "last_iteration": w.last_iteration, "temperature_pair": pair, "accepted": accepted,
                     "attempted": attempted, "rate": _rate(accepted, attempted)})
    rows.sort(key=lambda r: (r["window"], r["temperature_pair"]))
    return rows


def events(store, types=None):
    """Run events over both Phases, of `types` if given: iteration, phase, window, wall_time, type and payload (the
    JSON object parsed). Ordered by Sampling iteration, then Window, then as written."""
    relation = store.table("events", BOTH_PHASES)
    if relation is None:
        return []
    columns = ("iteration", "phase", "window", "wall_time", "type", "payload")
    rows = [dict(zip(columns, r)) for r in relation.select(", ".join(f'"{c}"' for c in columns)).fetchall()]
    rows = [{**r, "payload": json.loads(r["payload"])} for r in rows if types is None or r["type"] in types]
    rows.sort(key=lambda r: (r["iteration"], r["window"]))  # stable: events at one iteration keep their order
    return rows


def marked_events(store):
    """The events marked on the panel's time axes: non-finite values, interrupts, burn-in end, resumes and Run end."""
    return events(store, MARKED)


def ladder(store):
    """The temperature ladder at each `temperatures_adapted` event: iteration, phase, applied (false when the new
    ladder had NaNs and the old one was kept), temperatures (by temperature slot) and barriers (swap barrier per
    neighbouring pair). [] for a Run with one chain."""
    return [{"iteration": e["iteration"], "phase": e["phase"], **{k: e["payload"].get(k) for k in
                                                                ("applied", "temperatures", "barriers")}}
            for e in events(store, ("temperatures_adapted",))]


def burn_in_end(store):
    """The last `burn_in_ended` event: iteration, proposal_variances (sampler id to its final proposal variances, one
    per proposal scale; samplers without a scale are left out) and temperatures (the final ladder). None if burn-in
    has not ended in a closed Window."""
    ended = events(store, ("burn_in_ended",))
    if not ended:
        return None
    e = ended[-1]
    return {"iteration": e["iteration"], "proposal_variances": e["payload"].get("proposal_variances", {}),
            "temperatures": e["payload"].get("temperatures", [])}


def panel(store):
    """Everything the sampler panel draws for one Run: its closed Windows, acceptance, swaps, ladder, burn-in end and
    marked events."""
    windows = [{"window": w.number, "phase": w.phase, "first_iteration": w.first_iteration,
                "last_iteration": w.last_iteration} for w in store.closed_windows()]
    return {"windows": windows, "acceptance": acceptance(store), "swaps": swaps(store), "ladder": ladder(store),
            "burn_in_ended": burn_in_end(store), "events": marked_events(store)}
