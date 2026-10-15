"""pt-analysis: read Run stores through an Analysis manifest.

    pt-analysis info analysis.json      each Run's closed Windows per Phase and its Draw count
"""

import argparse
import sys

from pt_analysis.analysis import Analysis, AnalysisError
from pt_analysis.store import BOTH_PHASES


def describe_windows(windows):
    """`2 Windows (1-2), iterations 21-70` for a Phase's closed Windows, or `0 Windows`."""
    if not windows:
        return "0 Windows"
    numbers = [w.number for w in windows]
    span = str(numbers[0]) if len(numbers) == 1 else f"{numbers[0]}-{numbers[-1]}"
    noun = "Window" if len(windows) == 1 else "Windows"
    return f"{len(windows)} {noun} ({span}), iterations {windows[0].first_iteration}-{windows[-1].last_iteration}"


def info(analysis):
    print(f"Analysis: {analysis.path}")
    print(f"input: {analysis.input}")
    print(f"Runs: {len(analysis.runs)}")
    for label, store in analysis.runs.items():
        print(f"\n{label}: {store.path}")
        for phase in BOTH_PHASES:
            print(f"  {phase}: {describe_windows(store.windows((phase,)))}")
        print(f"  Draws: {store.draws()}")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="pt-analysis", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("info", help="list each Run's closed Windows per Phase and its Draws").add_argument("manifest")
    args = parser.parse_args(argv)
    try:
        analysis = Analysis.open(args.manifest)
    except AnalysisError as e:
        print(e, file=sys.stderr)
        return 1
    {"info": info}[args.command](analysis)
    return 0


if __name__ == "__main__":
    sys.exit(main())
