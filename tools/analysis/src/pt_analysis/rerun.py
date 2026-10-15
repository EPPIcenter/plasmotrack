"""Analysis manifests for the thesis rerun's layout (U7).

The rerun writes each fit's Chains to `<root>/<model>/<sim>/chain<N>`, one Run store each, where `<model>` is `full`
(every Infection), `partial` (the partially observed subset) or `null` (the full input under `--null-model`) and
`<sim>` is the simulation's path under the inputs directory, for example
`total_nodes_200_founder_prop_0.05_super_infected_prop_0_r0_1.5_mean_moi_2_asymptomatic_prop_0.5/simple_tx_network_med_div/sim_0_s7`.
Beside the engine inputs `<inputs>/<sim>_to_run.json.gz` and `<sim>_to_run_partially_obs.json.gz` is the
simulation's truth, `<inputs>/<sim>.json.gz`.

`manifests_from_rerun` writes `<out>/<model>/<sim>.json` per (model, simulation) with a Chain directory holding a
store (a `manifest.json`): the Chains in number order, the input, the IDPs, the truth, and the fit-summary factors
parsed from `<sim>` as the thesis pipeline's `parse_run_ids.R` does. The full fit names the null fit of the same
simulation as its null-model Analysis. Paths are relative to each manifest.
"""

import json
import os
import re
from pathlib import Path

INPUT_SUFFIX = {"full": "_to_run.json.gz", "partial": "_to_run_partially_obs.json.gz", "null": "_to_run.json.gz"}
# (missingness_condition, fit_type) of each model, as the per-fit summary names them.
FIT = {"full": ("full", "posterior"), "partial": ("missing_50", "posterior"), "null": ("full", "null_prior")}
PANEL_ALLELES = {"low": 5.0, "med": 10.0, "high": 20.0}
NUMBER = r"([0-9.]+)"
# The scenario directory, then the network directory: `[simple_]tx_network_<panel>_div` in the rerun's inputs, or
# `<panel>` in the pipeline's own layout, where `simple_` prefixes the scenario instead.
PATTERN = re.compile(
    rf"(?:^|/)(simple_)?total_nodes_{NUMBER}_founder_prop_{NUMBER}_super_infected_prop_{NUMBER}_r0_{NUMBER}"
    rf"_mean_moi_{NUMBER}_asymptomatic_prop_{NUMBER}/(simple_)?(?:tx_network_)?(low|med|high)(?:_div)?"
    r"/sim_([0-9]+)_s([0-9]+)$")


class RerunError(Exception):
    """A path does not follow the rerun's layout."""


def r_number(x):
    """A number as R's `paste0` writes it: integers without a decimal point."""
    return str(int(x)) if float(x).is_integer() else repr(float(x))


def parse_simulation(sim):
    """The simulation factors of a `<sim>` path, with its `run_key` (`parse_run_ids.R`)."""
    match = PATTERN.search(str(sim))
    if not match:
        raise RerunError(f"cannot read simulation factors from the path {sim}")
    simple, nodes, fp, sp, r0, coi, asym, simple_network, panel, index, seed = match.groups()
    model = "simplified" if simple or simple_network else "biological"
    factors = {"simulation_model": model, "total_nodes": float(nodes), "founder_probability": float(fp),
               "superinfection_probability": float(sp), "secondary_case_rate": float(r0),
               "founder_coi_rate": float(coi), "asymptomatic_probability": float(asym), "panel_name": panel,
               "panel_alleles_per_locus": PANEL_ALLELES[panel], "simulation_index": int(index), "seed": int(seed)}
    run_key = "__".join([model, f"fp{r_number(fp)}", f"sp{r_number(sp)}", f"r0{r_number(r0)}", f"coi{r_number(coi)}",
                         panel, f"sim{index}", f"s{seed}"])
    return {"run_key": run_key, **factors}


def chain_number(path):
    return int(path.name.removeprefix("chain"))


def fits(root):
    """{(model, sim): [Chain store directory, ...]} under `root`, Chains in number order; directories without a store
    are left out, and so is a fit with none."""
    found = {}
    for model in INPUT_SUFFIX:
        base = root / model
        if not base.is_dir():
            continue
        for chain in base.rglob("chain*"):
            if chain.is_dir() and re.fullmatch(r"chain[0-9]+", chain.name) and (chain / "manifest.json").is_file():
                found.setdefault((model, chain.parent.relative_to(base).as_posix()), []).append(chain)
    return {key: sorted(chains, key=chain_number) for key, chains in sorted(found.items())}


def manifests_from_rerun(root, inputs, out, symptomatic_idp=None, asymptomatic_idp=None):
    """Write one Analysis manifest per (model, simulation) fit under `root`; return their paths. IDPs default to
    `symptomatic_idp.csv` and `asymptomatic_idp.csv` beside the inputs directory, as the rerun keeps them."""
    root, inputs, out = Path(root).resolve(), Path(inputs).resolve(), Path(out).resolve()
    symptomatic_idp = Path(symptomatic_idp or inputs.parent / "symptomatic_idp.csv").resolve()
    asymptomatic_idp = Path(asymptomatic_idp or inputs.parent / "asymptomatic_idp.csv").resolve()
    found = fits(root)
    written = []
    for (model, sim), chains in found.items():
        path = out / model / f"{sim}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        rel = lambda p: os.path.relpath(p, path.parent)  # noqa: E731
        factors = parse_simulation(sim)
        missingness, fit_type = FIT[model]
        factors |= {"missingness_condition": missingness, "fit_type": fit_type}
        manifest = {"label": f"{factors['run_key']} {model}", "input": rel(inputs / f"{sim}{INPUT_SUFFIX[model]}"),
                    "symptomatic_idp": rel(symptomatic_idp), "asymptomatic_idp": rel(asymptomatic_idp),
                    "runs": {c.name: rel(c) for c in chains}, "truth": rel(inputs / f"{sim}.json.gz"),
                    "factors": factors}
        if model == "full" and ("null", sim) in found:
            manifest["null_model"] = rel(out / "null" / f"{sim}.json")
        path.write_text(json.dumps(manifest, indent=2) + "\n")
        written.append(path)
    return written
