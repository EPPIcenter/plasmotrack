"""The sampler panel's data (U4), for one Run: are the samplers doing their jobs?

Blocks (`acceptance`, `swaps`) hold counts per Window, not per Draw, so they are read over both Phases, and each row
carries its Window's Phase and Sampling-iteration range. Run events come from `events`, whose payloads are JSON
objects (`doc/store-layout.md`).
"""

import json
from statistics import median

from pt_analysis.sql import Scope
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


def _columns(names, rows):
    """`rows` as one list per field: {name: [value, ...]}, the payload shape for many rows."""
    return {name: list(values) for name, values in zip(names, zip(*rows))} if rows else {name: [] for name in names}


def _typed(s, store):
    """SQL over the Run's acceptance rows (both Phases) with each sampler's type and target and its rate (null when
    nothing was proposed), a sampler missing from `structure/samplers` split at its first `/`."""
    return f"""
        SELECT a."window", a.sampler, a.accepted, a.proposed,
               coalesce(k."type", split_part(a.sampler, '/', 1)) AS "type",
               coalesce(k.target, CASE WHEN contains(a.sampler, '/')
                                       THEN a.sampler[strpos(a.sampler, '/') + 1:] ELSE '' END) AS target,
               CASE WHEN a.proposed > 0 THEN a.accepted / a.proposed END AS rate
        FROM {s.name(store.table("acceptance", BOTH_PHASES))} a
             LEFT JOIN {s.name(store.structure("samplers"))} k USING (sampler)"""


BY_TYPE = ("type", "window", "targets", "accepted", "proposed", "rate", "rate_min", "rate_median", "rate_max")


def acceptance_by_type(store):
    """Acceptance per sampler type per Window, columnar ({field: [value per row]}): type, window, targets (samplers
    of the type with a row in the Window), accepted and proposed summed over them, rate (their ratio; None when
    nothing was proposed), and rate_min, rate_median and rate_max over the targets' own rates (targets that proposed
    nothing have none; None when no target has one). Ordered by type, then Window."""
    if store.table("acceptance", BOTH_PHASES) is None:
        return _columns(BY_TYPE, [])
    with Scope(store.connection) as s:
        rows = s.rows(f"""
            SELECT "type", "window", count(*), sum(accepted)::BIGINT, sum(proposed)::BIGINT,
                   CASE WHEN sum(proposed) > 0 THEN sum(accepted) / sum(proposed) END,
                   min(rate), median(rate), max(rate)
            FROM ({_typed(s, store)}) GROUP BY ALL ORDER BY "type", "window"
            """)
    return _columns(BY_TYPE, rows)


LIMIT = 500
PER_TARGET = ("target", "sampler", "window", "accepted", "proposed", "rate")


def acceptance_targets(store, kind, limit=LIMIT):
    """One sampler type's acceptance per target per Window, its problem targets first: type, targets (how many the
    type has), shown (how many are given: at most `limit`) and acceptance, columnar: target, sampler, window,
    accepted, proposed and rate (None when nothing was proposed). Targets are ordered by their mean rate over
    Windows, lowest first (targets that never proposed anything last), then by target; rows by target, then Window."""
    if store.table("acceptance", BOTH_PHASES) is None:
        return {"type": kind, "targets": 0, "shown": 0, "acceptance": _columns(PER_TARGET, [])}
    with Scope(store.connection) as s:
        rows = s.rows(f"""
            WITH typed AS (SELECT * FROM ({_typed(s, store)}) WHERE "type" = $kind),
                 ranked AS (SELECT target, row_number() OVER (ORDER BY avg(rate) ASC NULLS LAST, target) AS place
                            FROM typed GROUP BY target)
            SELECT t.target, t.sampler, t."window", t.accepted, t.proposed, t.rate, (SELECT count(*) FROM ranked)
            FROM typed t JOIN ranked r USING (target)
            WHERE r.place <= $limit ORDER BY r.place, t."window"
            """, {"kind": kind, "limit": limit})
        total = rows[0][-1] if rows else 0
    columns = _columns(PER_TARGET, [r[:-1] for r in rows])
    return {"type": kind, "targets": total, "shown": len(set(columns["target"])), "acceptance": columns}


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


def variance_summary(variances):
    """Burn-in-end proposal variances (sampler id to its variances) summarised per sampler type, columnar: type,
    samplers, scales (variances over them), and the min, median and max variance. Ordered by type."""
    by_type = {}
    for name, values in variances.items():
        entry = by_type.setdefault(_split(name)[0], [0, []])
        entry[0] += 1
        entry[1].extend(values)
    rows = [(kind, n, len(v), min(v), median(v), max(v)) for kind, (n, v) in sorted(by_type.items()) if v]
    return _columns(("type", "samplers", "scales", "min", "median", "max"), rows)


def panel(store):
    """Everything the sampler panel draws for one Run: its closed Windows, acceptance per sampler type
    (`acceptance_by_type`), swaps, ladder, burn-in end with its proposal variances per type (`variance_summary`) and
    marked events. The burn_in_ended event's payload leaves out its per-sampler proposal variances (one entry per
    sampler: `burn_in_ended` has them per type)."""
    windows = [{"window": w.number, "phase": w.phase, "first_iteration": w.first_iteration,
                "last_iteration": w.last_iteration} for w in store.closed_windows()]
    end = burn_in_end(store)
    if end is not None:
        end = {**end, "proposal_variances": variance_summary(end["proposal_variances"])}
    return {"windows": windows, "acceptance": acceptance_by_type(store), "swaps": swaps(store),
            "ladder": ladder(store), "burn_in_ended": end,
            "events": [{**e, "payload": {k: v for k, v in e["payload"].items() if k != "proposal_variances"}}
                       if e["type"] == "burn_in_ended" else e for e in marked_events(store)]}


def type_panel(store, kind, limit=LIMIT):
    """One sampler type in detail: its targets' acceptance (`acceptance_targets`, at most `limit` targets, problem
    targets first) and the burn-in-end proposal variances of the targets shown, in their order, columnar: target,
    sampler and variances (one per proposal scale). Targets without a proposal scale are left out."""
    targets = acceptance_targets(store, kind, limit)
    end = burn_in_end(store)
    variances = end["proposal_variances"] if end else {}
    order = dict.fromkeys(targets["acceptance"]["sampler"])
    rows = [(_split(name)[1], name, variances[name]) for name in order if name in variances]
    return {**targets, "proposal_variances": _columns(("target", "sampler", "variances"), rows)}
