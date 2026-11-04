#!/usr/bin/env python3
"""Convergence summary of independent chains written to RUN_DIR/chain*/ by the model CLI.

    chain_summary.py RUN_DIR [--chains 1,2,3] [--pairs A:B ...] [--show 15]

Treats the Run stores RUN_DIR/chain<K> (all of them, or those numbered in --chains) as the Chains of one Analysis,
writes its manifest to RUN_DIR/chain_summary.json and prints `pt-analysis summary` of it: per group of quantities,
the multi-chain effective sample size and split R-hat; for each pair A:B, how often each Chain had A infected before B
and how often that order changed between Draws; and, per Infection, its most probable parent pooled over Chains and
the largest difference between Chains in any of its parent probabilities. The derived tables are written under each
Chain's derived/ and RUN_DIR/chain_summary.derived/. The input and IDP paths in the manifest are those the first
Chain's Run manifest records. Needs the tools/analysis environment:
`uv run --project tools/analysis python tools/benchmark/chain_summary.py ...`.
"""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis" / "src"))
from pt_analysis.cli import main as pt_analysis  # noqa: E402

FLAGS = {"-i": "input", "--input": "input", "--symptomatic-idp": "symptomatic_idp",
         "--asymptomatic-idp": "asymptomatic_idp"}


def engine_inputs(run):
    """The input and IDP paths in a Run manifest's command line: `-i`/`--input` or the first argument, and the IDPs."""
    args = json.loads((run / "manifest.json").read_text()).get("cli_arguments", [])[1:]
    found = {"input": args[0]} if args and not args[0].startswith("-") else {}
    found |= {FLAGS[flag]: value for flag, value in zip(args, args[1:]) if flag in FLAGS}
    return {key: found.get(key, "") for key in ("input", "symptomatic_idp", "asymptomatic_idp")}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--chains", help="comma-separated chain numbers (default: all)")
    parser.add_argument("--pairs", nargs="*", default=[], help="infection pairs A:B")
    parser.add_argument("--show", type=int, default=15, help="infections listed, least agreement first")
    args = parser.parse_args()
    run = Path(args.run_dir)
    chains = sorted((p for p in run.iterdir() if p.is_dir() and re.fullmatch(r"chain\d+", p.name)),
                    key=lambda p: int(p.name[5:]))
    if args.chains:
        wanted = {f"chain{c}" for c in args.chains.split(",")}
        chains = [c for c in chains if c.name in wanted]
    if not chains:
        sys.exit(f"no chain<K> directories in {run}")
    manifest = run / "chain_summary.json"
    manifest.write_text(json.dumps({**engine_inputs(chains[0]), "runs": {c.name: c.name for c in chains}}, indent=2))
    return pt_analysis(["summary", str(manifest), "--show", str(args.show), "--pairs", *args.pairs])


if __name__ == "__main__":
    sys.exit(main())
