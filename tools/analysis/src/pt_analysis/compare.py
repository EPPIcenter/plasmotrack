"""Comparison of Analyses (U5): several Analyses of one input, opened together so their convergence can be read side
by side. The first is the primary Analysis, the one every single-Analysis view shows.

They are named on the command line (`pt-analysis serve a.json b.json`) or by a comparison manifest:

    {"analyses": ["build-a/analysis.json", "build-b/analysis.json"]}     paths resolved against the manifest's directory

Cost-normalised ESS is reported for an Analysis whose manifest carries `timings` in the shape of
`mixing_benchmark.py`'s `timings.json`, keyed by Chain label:
`{"chain1": {"sampling_seconds": 812.4, "instructions_retired": 3.1e12, ...}, ...}`. Each field is a cost when every
Chain has it: the sum over the Chains of `sampling_seconds` gives ESS per second (the pooled ESS over that sum), and
the sum of `instructions_retired` gives ESS per instruction retired, the measure the benchmark protocol uses (the
count is measured by the user, for example with `perf stat -e instructions`). An Analysis may have either or both.

Every Run of every Analysis must have been fitted to the same input: the `input` entry of `input_sha256` in each Run
manifest must agree, or the comparison is refused with `DifferentInputs`.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from pt_analysis.analysis import Analysis, AnalysisError
from pt_analysis.derive import GROUPS, derived_analysis


class DifferentInputs(AnalysisError):
    """The Analyses, or the Runs of one Analysis, were not fitted to the same input."""


@dataclass
class Comparison:
    analyses: dict[str, Analysis]  # label -> Analysis, primary first
    input_sha256: str | None  # None only for one Analysis whose Runs record no single input hash

    @property
    def primary(self):
        return next(iter(self.analyses.values()))


def manifest_paths(paths):
    """Analysis manifest paths from `paths`, each an Analysis manifest or a comparison manifest (one with
    `analyses`), expanded in place."""
    out = []
    for path in map(Path, paths):
        try:
            manifest = json.loads(path.read_text())
        except FileNotFoundError:
            raise AnalysisError(f"{path}: manifest not found") from None
        if isinstance(manifest, dict) and "analyses" in manifest:
            if not isinstance(manifest["analyses"], list) or not manifest["analyses"]:
                raise AnalysisError(f"{path}: 'analyses' must be a non-empty list of Analysis manifest paths")
            out.extend(path.resolve().parent / p for p in manifest["analyses"])
        else:
            out.append(path)
    return out


def input_hash(analysis):
    """The input hash every Run of `analysis` records; fails naming the Runs when they disagree or record none."""
    hashes = {label: store.manifest.get("input_sha256", {}).get("input") for label, store in analysis.runs.items()}
    if None in hashes.values():
        missing = ", ".join(label for label, h in hashes.items() if h is None)
        raise AnalysisError(f"{analysis.path}: Run manifest of {missing} records no input_sha256")
    if len(set(hashes.values())) > 1:
        raise DifferentInputs(f"{analysis.path}: its Runs were fitted to different inputs: "
                              + ", ".join(f"{label} {h[:12]}" for label, h in hashes.items()))
    return next(iter(hashes.values()))


def unique(label, taken):
    """`label`, or `label (2)`, `label (3)`, ... if it is taken."""
    candidate, n = label, 1
    while candidate in taken:
        n += 1
        candidate = f"{label} ({n})"
    return candidate


def open_comparison(paths):
    """Open the Analyses named by `paths` (see the module docstring), the first primary. Labels repeated between
    Analyses get a suffix ` (2)`, ` (3)`. Raises `DifferentInputs` naming each Analysis and its input hash when they
    were not all fitted to the same input; one Analysis alone is not checked."""
    analyses = {}
    for path in manifest_paths(paths):
        analysis = Analysis.open(path)
        analysis.label = unique(analysis.label, analyses)
        analyses[analysis.label] = analysis
    if len(analyses) == 1:  # nothing to compare: opened as `Analysis.open` opens it, hash or not
        try:
            return Comparison(analyses=analyses, input_sha256=input_hash(next(iter(analyses.values()))))
        except AnalysisError:
            return Comparison(analyses=analyses, input_sha256=None)
    hashes = {label: input_hash(a) for label, a in analyses.items()}
    if len(set(hashes.values())) > 1:
        raise DifferentInputs("Analyses of different inputs cannot be compared: "
                              + "; ".join(f"{label} ({analyses[label].path}) input {h[:12]}"
                                          for label, h in hashes.items()))
    return Comparison(analyses=analyses, input_sha256=next(iter(hashes.values())))


def chain_cost(analysis, field):
    """The sum over the Chains of `field` in the Analysis's timings; None without timings or when a Chain has no
    number in `field` there."""
    timings = analysis.timings
    if not isinstance(timings, dict):
        return None
    costs = [(timings.get(label) or {}).get(field) for label in analysis.runs]
    if any(isinstance(c, bool) or not isinstance(c, (int, float)) for c in costs):
        return None
    return float(sum(costs))


def sampling_seconds(analysis):
    """The sum over the Chains of `sampling_seconds` in the Analysis's timings (see `chain_cost`)."""
    return chain_cost(analysis, "sampling_seconds")


def instructions_retired(analysis):
    """The sum over the Chains of `instructions_retired` in the Analysis's timings (see `chain_cost`)."""
    return chain_cost(analysis, "instructions_retired")


# Cell key of cost-normalised ESS -> the cost it divides by.
PER_COST = {"ess_per_second": sampling_seconds, "ess_per_instruction": instructions_retired}


def compare_convergence(analyses):
    """The pooled rows of the convergence table of each of `analyses` (label -> Analysis, as `Comparison.analyses`),
    matched by quantity: one row per (group, name) in table order (groups in display order, names sorted),
    `{group, name, analyses: {label: {ess, rhat}}}`, every Analysis present and its ESS and R-hat null where it lacks
    the quantity. An Analysis with `sampling_seconds` timings also has `ess_per_second` in its cells, and one with
    `instructions_retired` timings `ess_per_instruction`."""
    costs = {label: {key: cost(a) for key, cost in PER_COST.items() if cost(a) is not None}
             for label, a in analyses.items()}
    cells = {}
    for label, analysis in analyses.items():
        rows = derived_analysis(analysis, "convergence").filter("chain IS NULL").fetchall()
        for group, name, _, ess, rhat in rows:
            cell = {"ess": ess, "rhat": rhat}
            for key, cost in costs[label].items():
                cell[key] = ess / cost if ess is not None and cost > 0 else None
            cells.setdefault((group, name), {})[label] = cell

    def empty(label):
        return {"ess": None, "rhat": None, **dict.fromkeys(costs[label])}

    order = sorted(cells, key=lambda key: (GROUPS.index(key[0]) if key[0] in GROUPS else len(GROUPS), key))
    return [{"group": group, "name": name,
             "analyses": {label: cells[(group, name)].get(label) or empty(label) for label in analyses}}
            for group, name in order]
