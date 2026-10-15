"""An Analysis: a set of Runs of one input, treated as Chains of one posterior, plus optional truth and epi data.

Described by a JSON manifest the user writes or a loader generates:

    {
      "input": "sim_0_s7_to_run.json.gz",            engine input the Runs were fitted to
      "symptomatic_idp": "symptomatic_idp.txt",      IDP files given to the engine
      "asymptomatic_idp": "asymptomatic_idp.txt",
      "runs": {"chain1": "chain1", "chain2": "chain2"},   label -> Run store directory, in display order
      "truth": "sim_0_s7.json.gz",                   optional: the simulation's true network, genotypes and times
      "epi": "epi.csv",                              optional: CSV keyed by Infection id, one row per Infection
      "pairwise_scores": "ibd.csv",                  optional: pairwise scores CSV (IBD), for the ranking baseline
      "null_model": "../null/analysis.json",         optional: the null-model Analysis of the same input
      "timings": "timings.json"                      optional: a timings JSON file, or the object inline
    }

Paths are resolved against the manifest's directory. Every Run directory must exist; the other files are named, not
checked, so an Analysis can be described before its truth or epi data is on disk.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

from pt_analysis.store import Store


class AnalysisError(Exception):
    """The manifest cannot be opened as an Analysis."""


@dataclass
class Analysis:
    path: Path
    input: Path
    symptomatic_idp: Path
    asymptomatic_idp: Path
    runs: dict[str, Store] = field(default_factory=dict)
    truth: Path | None = None
    epi: Path | None = None
    pairwise_scores: Path | None = None
    null_model: Path | None = None
    timings: dict | None = None

    @classmethod
    def open(cls, path):
        """Open the manifest at `path`. Fails naming the path of any Run directory that does not exist."""
        path = Path(path).resolve()
        manifest = json.loads(path.read_text())
        base = path.parent

        def resolve(key, required=False):
            if key not in manifest:
                if required:
                    raise AnalysisError(f"{path}: manifest has no '{key}'")
                return None
            return (base / manifest[key]).resolve()

        if not isinstance(manifest.get("runs"), dict) or not manifest["runs"]:
            raise AnalysisError(f"{path}: manifest needs a non-empty 'runs' object of label -> Run directory")
        runs = {}
        for label, run in manifest["runs"].items():
            run_dir = (base / run).resolve()
            if not run_dir.is_dir():
                raise AnalysisError(f"{path}: Run '{label}' directory not found: {run_dir}")
            runs[label] = Store(run_dir)

        timings = manifest.get("timings")
        if isinstance(timings, str):
            timings = json.loads((base / timings).read_text())

        return cls(path=path, input=resolve("input", required=True),
                   symptomatic_idp=resolve("symptomatic_idp", required=True),
                   asymptomatic_idp=resolve("asymptomatic_idp", required=True), runs=runs,
                   truth=resolve("truth"), epi=resolve("epi"), pairwise_scores=resolve("pairwise_scores"),
                   null_model=resolve("null_model"), timings=timings)
