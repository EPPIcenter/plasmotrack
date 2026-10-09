#!/usr/bin/env python3
"""Mixing benchmark: effective samples per second of sampling, from a burned-in state.

    mixing_benchmark.py run --model BIN --input JSON --symptomatic-idp F --asymptomatic-idp F --start DIR --out DIR
                            [--runs K --burnin B --samples S --seed N]
    mixing_benchmark.py analyze --out DIR

`run` copies the burned-in output directory `--start` once per run and hotloads each copy, so every run starts from
the same state. Hotloading does not restore the samplers' adapted proposal scales, so each run burns in for `--burnin`
sweeps (adaptation runs from sweep 20 until burn-in ends) before sampling `--samples` sweeps. Runs execute in
parallel, one per process; sweep times are taken from the model's progress lines as they are printed, so they are
wall-clock times and only comparable between benchmarks run on an otherwise idle machine.

`analyze` reads each run's traces and reports, per group of quantities, the multi-chain effective sample size across
runs (which counts disagreement between runs as variance), split R-hat and the share of series above 1.1, and
effective samples per second of sampling time summed over runs. Standard library only.
"""

import argparse
import gzip
import json
import math
import os
import pty
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path

PROGRESS = re.compile(r"^\((b|s)=(\d+)\)")


def run_one(args, run_dir, seed, timings):
    if run_dir.exists():
        shutil.rmtree(run_dir)
    # The start directory holds a file per infection, locus and parameter; on macOS clone it (copy-on-write) rather
    # than copying each file, which takes minutes.
    cloned = sys.platform == "darwin" and subprocess.run(["cp", "-cR", args.start, str(run_dir)]).returncode == 0
    if not cloned:
        shutil.copytree(args.start, run_dir)
    cmd = [args.model, "-i", args.input, "-o", str(run_dir), "-h",
           "--symptomatic-idp", args.symptomatic_idp, "--asymptomatic-idp", args.asymptomatic_idp,
           "-b", str(args.burnin), "-s", str(args.samples), "-t", "1", "--seed", str(seed)]
    # The model writes progress through stdio, which buffers whole blocks into a pipe; a pseudo-terminal makes it
    # flush each line, so lines can be timed as they arrive.
    stamps = []
    master, slave = pty.openpty()
    proc = subprocess.Popen(cmd, stdout=slave, stderr=slave, close_fds=True)
    os.close(slave)
    with open(run_dir / "run.log", "w", buffering=1) as log, os.fdopen(master, "r", errors="replace") as stream:
        try:
            for line in stream:
                log.write(line)
                match = PROGRESS.match(line)
                if match:
                    stamps.append((match.group(1), int(match.group(2)), time.time()))
        except OSError:
            pass  # the terminal closes when the model exits
    proc.wait()
    last_burnin = max((t for phase, _, t in stamps if phase == "b"), default=None)
    sample_times = [t for phase, _, t in stamps if phase == "s"]
    timings[run_dir.name] = {
        "exit": proc.returncode,
        "seed": seed,
        "sampling_seconds": (sample_times[-1] - last_burnin) if sample_times and last_burnin else None,
        "samples": len(sample_times),
    }


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    timings = {}
    threads = []
    for r in range(args.runs):
        thread = threading.Thread(target=run_one, args=(args, out / f"run{r}", args.seed + r, timings))
        thread.start()
        threads.append(thread)
    for thread in threads:
        thread.join()
    (out / "timings.json").write_text(json.dumps(timings, indent=2))
    print(json.dumps(timings, indent=2))


# --- traces ---------------------------------------------------------------------------------------------------------

def read_rows(path):
    with gzip.open(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def scalar_trace(path, skip_first):
    rows = read_rows(path)
    values = []
    for row in rows:
        try:
            values.append(float(row))
        except ValueError:
            continue  # header
    return values[1:] if skip_first else values


def vector_trace(path, skip_first):
    rows = [r for r in read_rows(path) if r[0].isdigit() or r[0] in "-."]
    values = [[float(x) for x in row.split(",")] for row in rows]
    return values[1:] if skip_first else values


def parent_set_traces(path, child):
    """Per iteration: the probability of the latent parent alone, and each observed parent's marginal probability."""
    latent = f"{{{child}_copy}}"
    no_parent = defaultdict(float)
    marginals = defaultdict(lambda: defaultdict(float))
    for row in read_rows(path)[1:]:
        members, prob, it = row.rsplit(",", 2)
        it, prob = int(it), float(prob)
        if members == latent:
            no_parent[it] += prob
        for member in members.strip("{}").split(";"):
            if not member.endswith("_copy"):
                marginals[member][it] += prob
    iters = sorted(set(no_parent) | {i for m in marginals.values() for i in m})
    return iters, [no_parent[i] for i in iters], {m: [v[i] for i in iters] for m, v in marginals.items()}


def run_traces(run_dir):
    """Every traced quantity of one run, as {group: {name: series}}. State traces start with the hotloaded state."""
    params = run_dir / "parameters"
    groups = defaultdict(dict)
    groups["log posterior"]["log posterior"] = scalar_trace(run_dir / "stats" / "likelihood.csv.gz", skip_first=False)
    for name in ["mean_coi", "mean_strains_tx", "parent_set_size_prob"]:
        groups["global scalars"][name] = scalar_trace(params / f"{name}.csv.gz", skip_first=True)
    for folder, group in [("eps_pos", "false positive rate"), ("eps_neg", "false negative rate"), ("infection_duration", "infection duration")]:
        for f in sorted((params / folder).glob("*.csv.gz")):
            groups[group][f.name] = scalar_trace(f, skip_first=True)
    for f in sorted((params / "allele_frequencies").glob("*.csv.gz")):
        for j, series in enumerate(zip(*vector_trace(f, skip_first=True))):
            groups["allele frequency"][f"{f.name}[{j}]"] = list(series)
    for f in sorted((run_dir / "parent_sets").glob("*_ps.csv.gz")):
        child = f.name[: -len("_ps.csv.gz")]
        _, no_parent, marginals = parent_set_traces(f, child)
        groups["P(no sampled parent)"][child] = no_parent
        for parent, series in marginals.items():
            groups["_edges"][f"{child}<-{parent}"] = series
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
        sys.exit(f"runs without a complete sampling phase: {', '.join(sorted(failed))}")
    seconds = sum(t["sampling_seconds"] for t in timings.values())
    draws = sum(t["samples"] for t in timings.values())

    # Edges: keep each child's most probable parent (by mean marginal over all runs).
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
    r.add_argument("--start", required=True, help="burned-in output directory to hotload")
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
