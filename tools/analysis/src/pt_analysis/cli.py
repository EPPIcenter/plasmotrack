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
    pt-analysis summarize analysis.json [--out row.parquet] [--edges [PATH]]
                                                         the fit-summary row of an Analysis (the 70 columns of
                                                         simulation_fit_summary.rds), and its per-edge rows
    pt-analysis cohort ROWS_DIR|GLOB ... --out cohort.parquet
                                                         concatenate fit-summary rows, and their per-edge rows
    pt-analysis manifests-from-rerun ROOT --inputs DIR --out DIR [--symptomatic-idp F] [--asymptomatic-idp F]
                                                         one Analysis manifest per fit of the thesis-rerun layout
                                                         ROOT/<model>/<sim>/chainN
    pt-analysis serve analysis.json [--port 8765] [--no-browser]
                                                         derive what is not current, then serve the analysis GUI
                                                         on 127.0.0.1 and open it in the browser

The summaries derive a Run's tables first when they are not current.
"""

import argparse
import sys
from pathlib import Path

from pt_analysis import cohort, fit_summary, functionals, rerun, summary
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


def summarize_command(args):
    analysis = Analysis.open(args.manifest)
    row = fit_summary.summarize(analysis)
    if args.out is None:
        if args.edges is not None:
            raise CommandError("--edges needs --out, or a path")
        table(["column", "value"], [[k, "" if v is None else v] for k, v in row.items()])
        return
    out = Path(args.out)
    fit_summary.write_row(row, out, analysis.path)
    print(f"{out}: the fit summary of {analysis.label}")
    if args.edges is not None:
        edges = Path(args.edges) if args.edges else cohort.edges_path(out)
        n = fit_summary.write_edges(analysis, edges)
        print(f"{edges}: {n} per-edge rows")


def cohort_command(args):
    n = cohort.concatenate(args.rows, args.out)
    edges = cohort.edges_path(args.out)
    print(f"{args.out}: {n} fit-summary rows" + (f"; per-edge rows in {edges}" if edges.is_file() else ""))


def rerun_command(args):
    written = rerun.manifests_from_rerun(args.root, args.inputs, args.out, args.symptomatic_idp, args.asymptomatic_idp)
    print(f"{args.out}: {len(written)} manifests")


def serve_command(args):
    from pt_analysis.server import serve

    serve(Analysis.open(args.manifest), port=args.port, open_browser=not args.no_browser)


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
    p = sub.add_parser("summarize", help="write the fit-summary row of an Analysis")
    p.add_argument("manifest")
    p.add_argument("--out", help="Parquet file for the row (default: print it)")
    p.add_argument("--edges", nargs="?", const="", default=None,
                   help="also write the per-edge rows (default path: the row file's stem + .edges.parquet)")
    p.set_defaults(handler=summarize_command)
    p = sub.add_parser("cohort", help="concatenate fit-summary rows into one Parquet file")
    p.add_argument("rows", nargs="+", help="directories of row files, row files, or globs")
    p.add_argument("--out", required=True, help="the cohort Parquet file")
    p.set_defaults(handler=cohort_command)
    p = sub.add_parser("manifests-from-rerun", help="write an Analysis manifest per fit of the thesis rerun")
    p.add_argument("root", help="the rerun's results: ROOT/<model>/<sim>/chainN")
    p.add_argument("--inputs", required=True, help="the inputs directory: <sim>_to_run.json.gz and <sim>.json.gz")
    p.add_argument("--out", required=True, help="directory for the manifests, <out>/<model>/<sim>.json")
    p.add_argument("--symptomatic-idp", help="default: symptomatic_idp.csv beside the inputs directory")
    p.add_argument("--asymptomatic-idp", help="default: asymptomatic_idp.csv beside the inputs directory")
    p.set_defaults(handler=rerun_command)
    p = sub.add_parser("serve", help="serve the analysis GUI over an Analysis")
    p.add_argument("manifest")
    p.add_argument("--port", type=int, default=8765, help="port on 127.0.0.1 (default 8765; 0 picks a free one)")
    p.add_argument("--no-browser", action="store_true", help="do not open the browser")
    p.set_defaults(handler=serve_command)
    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except (AnalysisError, CommandError, summary.SummaryError, cohort.CohortError, rerun.RerunError) as e:
        print(e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
