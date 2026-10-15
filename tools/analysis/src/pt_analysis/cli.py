"""pt-analysis: read Run stores through an Analysis manifest.

    pt-analysis info analysis.json                       each Run's closed Windows per Phase and its Draw count
    pt-analysis derive RUN_DIR | analysis.json           write the derived tables under each Run's derived/, and an
                                                         Analysis's convergence table under <manifest stem>.derived/
    pt-analysis edges RUN_DIR [--top N] [--unobserved]   the top pooled Edge marginals
    pt-analysis map RUN_DIR [--level 0.8]                the MAP graph, psi_MAP and the parent-set size range
    pt-analysis credible RUN_DIR INFECTION [--level 0.8] an Infection's credible collection of Parent sets
    pt-analysis summary analysis.json [--pairs A:B ...] [--show N]
                                                         the chain summary: ESS and R-hat per group, the order of
                                                         Infection pairs per Chain, and agreement on parents

The summaries derive a Run's tables first when they are not current.
"""

import argparse
import sys
from pathlib import Path

from pt_analysis import functionals, summary
from pt_analysis.analysis import Analysis, AnalysisError
from pt_analysis.derive import analysis_dir, derive, derive_analysis
from pt_analysis.store import BOTH_PHASES, Store


class CommandError(Exception):
    """A command cannot run on what it was given; the message says why."""


def describe_windows(windows):
    """`2 Windows (1-2), iterations 21-70` for a Phase's closed Windows, or `0 Windows`."""
    if not windows:
        return "0 Windows"
    numbers = [w.number for w in windows]
    span = str(numbers[0]) if len(numbers) == 1 else f"{numbers[0]}-{numbers[-1]}"
    noun = "Window" if len(windows) == 1 else "Windows"
    return f"{len(windows)} {noun} ({span}), iterations {windows[0].first_iteration}-{windows[-1].last_iteration}"


def open_run(path):
    path = Path(path).resolve()
    if not (path / "manifest.json").is_file():
        raise CommandError(f"{path}: not a Run store (no manifest.json)")
    return Store(path)


def table(header, rows):
    """`rows` under `header`, columns left-aligned and separated by two spaces."""
    rows = [header] + [[str(v) for v in row] for row in rows]
    widths = [max(len(row[i]) for row in rows) for i in range(len(header))]
    for row in rows:
        print("  ".join(v.ljust(w) for v, w in zip(row, widths)).rstrip())


def number(x):
    return f"{x:.4g}"


def members(ranked_set):
    return "{" + ", ".join(ranked_set.members) + "}"


def info(args):
    analysis = Analysis.open(args.manifest)
    print(f"Analysis: {analysis.path}")
    print(f"input: {analysis.input}")
    print(f"Runs: {len(analysis.runs)}")
    for label, store in analysis.runs.items():
        print(f"\n{label}: {store.path}")
        for phase in BOTH_PHASES:
            print(f"  {phase}: {describe_windows(store.windows((phase,)))}")
        print(f"  Draws: {store.draws()}")


def derive_command(args):
    path = Path(args.path)
    analysis = Analysis.open(path) if path.is_file() else None
    stores = list(analysis.runs.values()) if analysis else [open_run(path)]
    for store in stores:
        if derive(store):
            windows = ", ".join(f"{phase} Windows {[w.number for w in store.windows((phase,))]}"
                                for phase in BOTH_PHASES)
            print(f"{store.path}: derived from {windows}")
        else:
            print(f"{store.path}: up to date")
    if analysis:
        if derive_analysis(analysis):
            print(f"{analysis_dir(analysis)}: derived convergence from Chains {', '.join(analysis.runs)}")
        else:
            print(f"{analysis_dir(analysis)}: up to date")


def edges(args):
    rows = functionals.top_edges(open_run(args.run), n=args.top, unobserved=args.unobserved)
    table(["parent", "child", "weight"], [[e.parent, e.child, number(e.weight)] for e in rows])


def map_command(args):
    store = open_run(args.run)
    graph = functionals.map_graph(store)
    table(["infection", "set", "weight", "members"],
          [[s.infection, s.set_id, number(s.weight), members(s)] for s in graph.values()])
    psi = functionals.psi_map(store)
    lo, hi = functionals.size_range(store, args.level)
    print(f"psi_MAP: {number(psi)}")
    print(f"parent-set size range at {args.level:g}: [{number(lo)}, {number(hi)}]")


def credible(args):
    collection = functionals.credible_collection(open_run(args.run), args.infection, args.level)
    if not collection:
        raise CommandError(f"Infection {args.infection} has no Parent sets with weight in {args.run}")
    noun = "set" if len(collection) == 1 else "sets"
    print(f"Infection {args.infection}: {len(collection)} {noun} with mass at least {args.level:g}")
    table(["rank", "set", "weight", "cumulative", "class", "members"],
          [[s.rank, s.set_id, number(s.weight), number(s.cumulative_mass), s.composition, members(s)]
           for s in collection])


def summary_command(args):
    pairs = []
    for pair in args.pairs:
        a, sep, b = pair.partition(":")
        if not sep or not a or not b:
            raise CommandError(f"--pairs takes Infection pairs A:B, not '{pair}'")
        pairs.append((a, b))
    summary.print_summary(Analysis.open(args.manifest), pairs=pairs, show=args.show)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="pt-analysis", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("info", help="list each Run's closed Windows per Phase and its Draws")
    p.add_argument("manifest")
    p.set_defaults(handler=info)
    p = sub.add_parser("derive", help="write the derived tables of a Run, or of every Run of an Analysis")
    p.add_argument("path", help="a Run store directory or an Analysis manifest")
    p.set_defaults(handler=derive_command)
    p = sub.add_parser("edges", help="print the top pooled Edge marginals of a Run")
    p.add_argument("run")
    p.add_argument("--top", type=int, help="print only the N largest")
    p.add_argument("--unobserved", action="store_true", help="include edges from the unobserved source")
    p.set_defaults(handler=edges)
    p = sub.add_parser("map", help="print a Run's MAP graph, psi_MAP and parent-set size range")
    p.add_argument("run")
    p.add_argument("--level", type=float, default=functionals.LEVEL, help="mass for the size range (default 0.8)")
    p.set_defaults(handler=map_command)
    p = sub.add_parser("credible", help="print an Infection's credible collection of Parent sets")
    p.add_argument("run")
    p.add_argument("infection")
    p.add_argument("--level", type=float, default=functionals.LEVEL, help="mass to reach (default 0.8)")
    p.set_defaults(handler=credible)
    p = sub.add_parser("summary", help="print the chain summary of an Analysis")
    p.add_argument("manifest")
    p.add_argument("--pairs", nargs="*", default=[], help="Infection pairs A:B")
    p.add_argument("--show", type=int, default=15, help="Infections listed, least agreement first (default 15)")
    p.set_defaults(handler=summary_command)
    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except (AnalysisError, CommandError, summary.SummaryError) as e:
        print(e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
