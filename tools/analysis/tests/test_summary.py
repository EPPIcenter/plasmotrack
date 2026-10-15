"""`pt-analysis summary`: the chain summary of an Analysis, from its derived tables and the Ordering."""

import re
import shutil

import duckdb
import pytest
from conftest import GOLDEN
from test_analysis import write_manifest

from pt_analysis.analysis import Analysis
from pt_analysis.cli import main
from pt_analysis.summary import UNOBSERVED_ONLY, agreement, pair_order, spread_quantiles


def test_pair_order_counts_draws_with_a_first_and_switches_between_adjacent_draws():
    assert pair_order([True, True, False, True]) == (3, 4, 2)
    assert pair_order([False] * 5) == (0, 5, 0)
    assert pair_order([]) == (0, 0, 0)


CHAINS = [
    {"a": {"b": 0.9, "c": 0.1, UNOBSERVED_ONLY: 0.0}, "b": {UNOBSERVED_ONLY: 1.0},
     "c": {"a": 0.7, UNOBSERVED_ONLY: 0.3}},
    {"a": {"b": 0.2, UNOBSERVED_ONLY: 0.5}, "b": {UNOBSERVED_ONLY: 1.0}, "c": {"a": 0.6, UNOBSERVED_ONLY: 0.4}},
]


def test_agreement_pools_parent_probabilities_over_chains_and_sorts_by_spread():
    """Infection a: b pooled (0.9 + 0.2) / 2 = 0.55 with spread 0.7; c is missing from chain 2, so 0 there (spread
    0.1); unobserved only pooled 0.25 with spread 0.5. Infection b has no observed parent: top '-'."""
    rows = agreement(CHAINS)
    assert [(r.infection, r.top) for r in rows] == [("a", "b"), ("c", "a"), ("b", "-")]
    a, c, b = rows
    assert (a.p_top, a.p_unobserved_only, a.spread) == pytest.approx((0.55, 0.25, 0.7))
    assert (c.p_top, c.p_unobserved_only, c.spread) == pytest.approx((0.65, 0.35, 0.1))
    assert (b.p_top, b.p_unobserved_only, b.spread) == (0, 1.0, 0)


def test_spread_quantiles_take_the_upper_median_and_the_90th_percentile_by_index():
    assert spread_quantiles(agreement(CHAINS)) == pytest.approx((0.1, 0.7, 0.7))


@pytest.fixture
def analysis_path(tmp_path):
    for name in ("single-chain", "replica-exchange"):
        shutil.copytree(GOLDEN / name, tmp_path / name)
    return write_manifest(tmp_path / "analysis.json", runs={"chain1": "single-chain", "chain2": "replica-exchange"})


def columns(line):
    return re.split(r"\s{2,}", line.strip())


# chain_summary.py's table on the two golden stores as chain1 and chain2, before the port; the edge row is
# mixing_benchmark's ESS and R-hat over all three of run_traces' edge series.
BASELINE = """\
log posterior                  1         3.6       3.6       100%       1.26
global scalars                 3        26.0      19.5         0%      0.938
false positive rate            3        26.0      21.3         0%      0.999
false negative rate            3        26.0       3.2        33%       1.11
infection duration             3        26.0      24.1         0%       1.01
allele frequency              20        10.6       1.4        65%       2.43
edge                           3        26.0      12.5        33%        1.1"""


def test_the_diagnostics_table_reproduces_the_chain_summary(analysis_path, capsys):
    assert main(["summary", str(analysis_path)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "2 Chains: chain1, chain2"
    assert out[1] == ""
    assert columns(out[2]) == ["group", "series", "ESS median", "ESS min", "R-hat>1.1", "R-hat max"]
    table = [columns(line) for line in out[3:12]]
    assert table[:7] == [columns(line) for line in BASELINE.splitlines()]
    assert [row[:2] for row in table[7:]] == [["unobserved-source inclusion", "1"], ["functionals", "3"]]
    assert out[12] == ""


def test_the_agreement_section_on_the_golden_stores(analysis_path, capsys):
    assert main(["summary", str(analysis_path), "--show", "2"]) == 0
    out = capsys.readouterr().out.splitlines()
    start = out.index(next(line for line in out if columns(line)[0] == "infection"))
    assert columns(out[start]) == ["infection", "top parent", "P(top)", "P(unobserved only)", "spread"]
    # Inf 3: edge 1->3 pooled (0.89138 + 0.77570) / 2; unobserved only (0.1084 + 0.2241) / 2; the spread is edge
    # 1->3's 0.89138 - 0.77570.
    assert columns(out[start + 1]) == ["3", "1", "0.834", "0.166", "0.116"]
    assert columns(out[start + 2])[:2] == ["2", "1"]
    assert len(out) == start + 4
    assert out[start + 3].startswith("largest between-Chain difference in a parent probability, per Infection: ")
    assert out[start + 3].endswith("max 0.116")


def test_pairs_count_each_chains_draws_with_a_first_and_its_switches(analysis_path, capsys):
    """Swap Infections 1 and 2 in chain1's Ordering at Draws 35 and 40: 8 of 10 Draws have 1 before 2, and the order
    changes twice."""
    window = analysis_path.parent / "single-chain" / "ordering" / "window-000001.parquet"
    swapped = window.with_name("swapped.parquet")
    duckdb.sql(f"""COPY (SELECT * REPLACE (CASE WHEN iteration IN (35, 40) AND infection IN ('1', '2')
                                               THEN 1 - rank ELSE rank END AS rank)
                        FROM '{window}') TO '{swapped}' (FORMAT parquet)""")
    swapped.replace(window)
    assert main(["summary", str(analysis_path), "--pairs", "1:2", "2:3"]) == 0
    out = capsys.readouterr().out.splitlines()
    start = out.index(next(line for line in out if line.startswith("chain1: ")))
    assert out[start - 1] == ""
    assert out[start:start + 2] == ["chain1: 1 before 2: 8/10 (changes 2) | 2 before 3: 10/10 (changes 0)",
                                    "chain2: 1 before 2: 10/10 (changes 0) | 2 before 3: 10/10 (changes 0)"]
    assert out[start + 2] == ""


def test_a_pair_naming_an_unknown_infection_fails(analysis_path, capsys):
    assert main(["summary", str(analysis_path), "--pairs", "1:99"]) == 1
    assert "99" in capsys.readouterr().err


def test_the_summary_reads_derived_tables_it_derives_first(analysis_path):
    main(["summary", str(analysis_path)])
    analysis = Analysis.open(analysis_path)
    assert (analysis_path.parent / "analysis.derived" / "convergence.parquet").is_file()
    assert all((store.path / "derived" / "edges.parquet").is_file() for store in analysis.runs.values())
