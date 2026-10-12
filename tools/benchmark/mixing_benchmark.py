#!/usr/bin/env python3
"""Mixing benchmark: effective samples per second of sampling, from a burned-in state.

    mixing_benchmark.py run --model BIN --input JSON --symptomatic-idp F --asymptomatic-idp F --start CHECKPOINT --out DIR
                            [--runs K --burnin B --samples S --seed N]
    mixing_benchmark.py analyze --out DIR

`run` starts every Run with `--start-from` on the same checkpoint (`--start`): a Run's `checkpoint-burn-in.cbor`, or a
Run's latest `checkpoint.cbor`; a Run directory stands for its `checkpoint.cbor`. So every Run starts from the same
State and adaptation state. Each Run burns in for `--burnin` Sampling iterations with its own seed (adaptation runs from
Sampling iteration 20 until burn-in ends), then samples `--samples` Sampling iterations with `-t 1`. Runs execute in
parallel, one per process, each writing its store to DIR/runK. Sampling time is the wall time the store's
`window_closed` events record for the sampling Phase, so it is only comparable between benchmarks run on an otherwise
idle machine.

`analyze` reads each Run's store and reports, per group of quantities, the multi-chain effective sample size across
Runs (which counts disagreement between Runs as variance), split R-hat and the share of series above 1.1, and effective
samples per second of sampling time summed over Runs. It reads only closed Windows, those with a marker under
`_closed/`, so a killed Run is read up to its last closed Window. The store layout is described in doc/store-layout.md.

Needs pyarrow (`pip install pyarrow`); otherwise standard library only.
"""

import argparse
import json
import math
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def checkpoint_path(start):
    start = Path(start)
    return start / "checkpoint.cbor" if start.is_dir() else start


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    start = checkpoint_path(args.start)
    if not start.is_file():
        sys.exit(f"no checkpoint at {start}")
    started = []
    for r in range(args.runs):
        run_dir = out / f"run{r}"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir()
        seed = args.seed + r
        cmd = [args.model, "-i", args.input, "-o", str(run_dir), "--start-from", str(start),
               "--symptomatic-idp", args.symptomatic_idp, "--asymptomatic-idp", args.asymptomatic_idp,
               "-b", str(args.burnin), "-s", str(args.samples), "-t", "1", "--seed", str(seed)]
        # The console goes beside the store, so that the store holds only what the engine writes.
        console = open(out / f"run{r}.console", "w")
        started.append((run_dir, seed, console, subprocess.Popen(cmd, stdout=console, stderr=subprocess.STDOUT)))
    timings = {}
    for run_dir, seed, console, proc in started:
        proc.wait()
        console.close()
        timings[run_dir.name] = {"exit": proc.returncode, "seed": seed, **sampling_time(run_dir)}
    (out / "timings.json").write_text(json.dumps(timings, indent=2))
    print(json.dumps(timings, indent=2))


# --- the store ------------------------------------------------------------------------------------------------------

def closed_windows(run_dir):
    """The Windows a reader may use: those whose commit marker exists. In-flight files are dot-prefixed."""
    closed = run_dir / "_closed"
    if not closed.is_dir():
        return []
    return sorted(json.loads(p.read_text())["window"] for p in closed.glob("window-*"))


def read_table(run_dir, table, windows):
    """One table over the given closed Windows, or None if it has no rows in them. A Window in which a table has no
    rows has no file for it."""
    files = [run_dir / table / f"window-{w:06d}.parquet" for w in windows]
    parts = [pq.read_table(f) for f in files if f.exists()]
    return pa.concat_tables(parts) if parts else None


def sampling_rows(table):
    """The sampling-Phase rows of a Trace or of the events table, in Sampling iteration order."""
    if table is None:
        return []
    table = table.filter(pc.equal(table["phase"], "sampling"))
    return table.take(pc.sort_indices(table, [("iteration", "ascending")])).to_pylist()


def sampling_time(run_dir):
    """Wall time and Sampling iterations of the sampling Phase, summed over the `window_closed` events."""
    closes = [json.loads(e["payload"]) for e in sampling_rows(read_table(run_dir, "events", closed_windows(run_dir)))
              if e["type"] == "window_closed"]
    samples = sum(c["sampling_iterations"] for c in closes)
    return {"sampling_seconds": sum(c["seconds"] for c in closes) if samples else None, "samples": samples}


# --- traces ---------------------------------------------------------------------------------------------------------

def parent_set_traces(run_dir, windows, groups):
    """Per Sampling iteration: the probability of each Infection's Latent parent alone, and each observed parent's
    marginal probability. Parent sets below the Run's floor are not in the store and count as 0."""
    latent_parent = {r["infection"]: r["latent_parent"]
                     for r in pq.read_table(run_dir / "structure" / "infections.parquet").to_pylist()}
    latent_parents = set(latent_parent.values())
    members = defaultdict(list)  # Parent set members by (Infection, set_id), from every closed Window of both Phases
    member_rows = read_table(run_dir, "parent_set_members", windows)
    for r in member_rows.to_pylist() if member_rows is not None else []:
        members[(r["infection"], r["set_id"])].append(r["member"])

    iterations = defaultdict(set)
    no_parent = defaultdict(lambda: defaultdict(float))
    marginals = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    for r in sampling_rows(read_table(run_dir, "parent_sets", windows)):
        child, it, prob = r["infection"], r["iteration"], r["prob"]
        iterations[child].add(it)
        parents = members[(child, r["set_id"])]
        if parents == [latent_parent[child]]:
            no_parent[child][it] += prob
        for parent in parents:
            if parent not in latent_parents:
                marginals[child][parent][it] += prob
    for child in sorted(iterations):
        iters = sorted(iterations[child])
        groups["P(no sampled parent)"][child] = [no_parent[child][i] for i in iters]
        for parent, series in marginals[child].items():
            groups["_edges"][f"{child}<-{parent}"] = [series[i] for i in iters]


def run_traces(run_dir):
    """Every traced quantity of one Run's sampling Phase, as {group: {name: series}}."""
    windows = closed_windows(run_dir)
    groups = defaultdict(dict)

    scalars = sampling_rows(read_table(run_dir, "run_scalars", windows))
    groups["log posterior"]["log posterior"] = [r["prior"] + r["likelihood"] for r in scalars]
    for name in ["mean_coi", "mean_strains_tx", "parent_set_size_prob"]:
        groups["global scalars"][name] = [r[name] for r in scalars]

    infections = sampling_rows(read_table(run_dir, "infections", windows))
    for column, group in [("expected_false_positives", "false positive rate"),
                          ("expected_false_negatives", "false negative rate"),
                          ("infection_duration", "infection duration")]:
        for r in infections:
            groups[group].setdefault(r["infection"], []).append(r[column])

    for r in sampling_rows(read_table(run_dir, "allele_frequencies", windows)):
        groups["allele frequency"].setdefault(f"{r['locus']}[{r['allele']}]", []).append(r["frequency"])

    parent_set_traces(run_dir, windows, groups)
    return groups


# --- diagnostics ----------------------------------------------------------------------------------------------------

def autocovariances(x):
    n = len(x)
    mean = sum(x) / n
    centred = [v - mean for v in x]
    return [sum(centred[i] * centred[i + lag] for i in range(n - lag)) / n for lag in range(n)]


def multichain_ess(chains):
    """Effective sample size of several chains of equal length, counting disagreement between chains as variance
    (Vehtari et al. 2021, without rank normalisation): autocorrelations are taken against the pooled variance and
    summed with Geyer's initial monotone sequence. None for a constant series."""
    m = len(chains)
    n = min(len(c) for c in chains)
    if n < 4:
        return None
    chains = [c[:n] for c in chains]
    covs = [autocovariances(c) for c in chains]
    means = [sum(c) / n for c in chains]
    within = statistics.mean(cov[0] * n / (n - 1) for cov in covs)
    between = n * statistics.variance(means) if m > 1 else 0.0
    pooled = (n - 1) / n * within + between / n
    if pooled <= 0:
        return None
    rho = [1.0 - (within - statistics.mean(cov[t] for cov in covs)) / pooled for t in range(n)]
    total = 0.0
    previous = float("inf")
    for k in range(0, n - 1, 2):
        pair = rho[k] + rho[k + 1]
        if pair <= 0:
            break
        pair = min(pair, previous)  # initial monotone sequence
        total += pair
        previous = pair
    tau = max(-1.0 + 2.0 * total, 1.0 / math.log10(max(m * n, 10)))
    return m * n / tau


def split_rhat(chains):
    halves = []
    for c in chains:
        h = len(c) // 2
        if h >= 2:
            halves += [c[:h], c[h: 2 * h]]
    if len(halves) < 4:
        return None
    n = len(halves[0])
    means = [sum(h) / n for h in halves]
    within = statistics.mean(statistics.variance(h) for h in halves)
    between = n * statistics.variance(means)
    if within <= 0:
        # Constant within every half: identical everywhere (no information), or stuck at different values.
        return None if between <= 0 else math.inf
    return math.sqrt(((n - 1) / n * within + between / n) / within)


def analyze(args):
    out = Path(args.out)
    timings = json.loads((out / "timings.json").read_text())
    run_dirs = sorted(p for p in out.glob("run*") if p.is_dir())
    traces = [run_traces(d) for d in run_dirs]
    failed = [name for name, t in timings.items() if t["exit"] != 0 or t["sampling_seconds"] is None]
    if failed:
        sys.exit(f"Runs without a complete sampling Phase: {', '.join(sorted(failed))}")
    seconds = sum(t["sampling_seconds"] for t in timings.values())
    draws = sum(t["samples"] for t in timings.values())

    # Edges: keep each child's most probable parent (by mean marginal over all Runs).
    best = {}
    for run_groups in traces:
        for key, series in run_groups.get("_edges", {}).items():
            child = key.split("<-")[0]
            mean = sum(series) / len(series)
            best.setdefault(child, defaultdict(float))[key] += mean
    top_edges = {max(scores, key=scores.get) for scores in best.values()}
    for run_groups in traces:
        run_groups["top parent edge"] = {k: v for k, v in run_groups.pop("_edges", {}).items() if k in top_edges}

    print(f"{len(run_dirs)} runs, {draws} draws, {seconds:.0f} s of sampling ({seconds / max(draws, 1):.1f} s per sweep)\n")
    print(f"{'group':<24}{'series':>8}{'ESS median':>12}{'ESS min':>10}{'ESS/s median':>14}{'ESS/s min':>11}{'R-hat>1.1':>11}{'R-hat max':>11}")
    report = {}
    for group in traces[0]:
        names = sorted(set().union(*(set(t.get(group, {})) for t in traces)))
        totals, rhats = [], []
        for name in names:
            chains = [t[group][name] for t in traces if name in t.get(group, {})]
            total = multichain_ess(chains)
            if total is None:
                continue
            totals.append(total)
            r = split_rhat(chains)
            if r is not None:
                rhats.append(r)
        if not totals:
            continue
        report[group] = {"series": len(totals), "ess_median": statistics.median(totals), "ess_min": min(totals),
                         "ess_per_s_median": statistics.median(totals) / seconds, "ess_per_s_min": min(totals) / seconds,
                         "rhat_max": max(rhats) if rhats else None,
                         "rhat_above_1_1": sum(r > 1.1 for r in rhats) / len(rhats) if rhats else None}
        g = report[group]
        rhat = f"{g['rhat_max']:.3g}" if g["rhat_max"] is not None else "-"
        above = f"{g['rhat_above_1_1']:.0%}" if g["rhat_above_1_1"] is not None else "-"
        print(f"{group:<24}{g['series']:>8}{g['ess_median']:>12.1f}{g['ess_min']:>10.1f}{g['ess_per_s_median']:>14.4f}{g['ess_per_s_min']:>11.4f}{above:>11}{rhat:>11}")
    (out / "report.json").write_text(json.dumps({"draws": draws, "sampling_seconds": seconds, "groups": report}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--model", required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--symptomatic-idp", required=True)
    r.add_argument("--asymptomatic-idp", required=True)
    r.add_argument("--start", required=True,
                   help="checkpoint every Run starts from: checkpoint-burn-in.cbor or checkpoint.cbor, or a Run directory for its checkpoint.cbor")
    r.add_argument("--out", required=True)
    r.add_argument("--runs", type=int, default=6)
    r.add_argument("--burnin", type=int, default=40)
    r.add_argument("--samples", type=int, default=80)
    r.add_argument("--seed", type=int, default=1000)
    a = sub.add_parser("analyze")
    a.add_argument("--out", required=True)
    args = parser.parse_args()
    {"run": run, "analyze": analyze}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
