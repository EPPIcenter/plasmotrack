"""`summarize` and `cohort` from the command line: one row per Analysis, rows concatenated with their per-edge rows,
and the cohort read back for the cohort view."""

import pyarrow.parquet as pq
import pytest
from test_fit_summary import FACTORS
from test_metrics import open_analysis, stores  # noqa: F401  (fixture)

from pt_analysis import cohort, fit_summary
from pt_analysis.cli import main


@pytest.fixture
def rows(tmp_path, stores):  # noqa: F811
    """Two fit-summary rows with their edge files, as `summarize` writes them, under rows/."""
    (tmp_path / "rows").mkdir()
    for fit_type in ("posterior", "null_prior"):
        manifest = open_analysis(tmp_path, stores, name=f"{fit_type}.json", factors=FACTORS | {"fit_type": fit_type})
        assert main(["summarize", str(manifest.path), "--out", str(tmp_path / "rows" / f"{fit_type}.parquet"),
                     "--edges"]) == 0
    return tmp_path / "rows"


def test_summarize_writes_one_row_and_its_edges(rows):
    table = pq.read_table(rows / "posterior.parquet")
    assert table.schema.names == fit_summary.SCHEMA.names and table.num_rows == 1
    edges = pq.read_table(rows / "posterior.edges.parquet")
    assert edges.schema == fit_summary.EDGE_SCHEMA and edges.num_rows == 3


def test_summarize_prints_the_row_without_out(tmp_path, stores, capsys):  # noqa: F811
    manifest = open_analysis(tmp_path, stores, factors=FACTORS)
    assert main(["summarize", str(manifest.path)]) == 0
    assert "estimated_average_outdegree" in capsys.readouterr().out


def test_cohort_concatenates_rows_and_edges_and_names_each_manifest(rows, tmp_path):
    out = tmp_path / "cohort" / "cohort.parquet"
    assert main(["cohort", str(rows), "--out", str(out)]) == 0
    table = pq.read_table(out)
    assert table.num_rows == 2
    assert table.schema.names == [*fit_summary.SCHEMA.names, "manifest"]
    assert sorted(table.column("fit_type").to_pylist()) == ["null_prior", "posterior"]
    # Manifest paths are relative to the cohort file, so the cohort moves with the fits.
    assert sorted(table.column("manifest").to_pylist()) == ["../null_prior.json", "../posterior.json"]
    edges = pq.read_table(cohort.edges_path(out))
    assert edges.num_rows == 6 and set(edges.column("fit_type").to_pylist()) == {"posterior", "null_prior"}


def test_cohort_takes_a_glob(rows, tmp_path):
    out = tmp_path / "c.parquet"
    assert main(["cohort", str(rows / "post*.parquet"), "--out", str(out)]) == 0
    assert pq.read_table(out).column("fit_type").to_pylist() == ["posterior"]


def test_cohort_refuses_a_file_that_is_not_a_fit_summary(rows, tmp_path, capsys):
    pq.write_table(pq.read_table(rows / "posterior.parquet").drop(["seed"]), rows / "broken.parquet")
    assert main(["cohort", str(rows), "--out", str(tmp_path / "c.parquet")]) == 1
    assert "broken.parquet" in capsys.readouterr().err


def test_cohort_with_no_rows_fails(tmp_path, capsys):
    assert main(["cohort", str(tmp_path / "nothing*.parquet"), "--out", str(tmp_path / "c.parquet")]) == 1
    assert "no fit-summary rows" in capsys.readouterr().err


def test_load_gives_columns_and_calibration_counts_per_fit(rows, tmp_path):
    out = tmp_path / "cohort.parquet"
    cohort.concatenate([rows], out)
    payload = cohort.load(out, bins=10)
    assert payload["n"] == 2
    names = [c["name"] for c in payload["columns"]]
    assert names[:3] == ["run_key", "simulation_model", "total_nodes"] and names[-1] == "manifest"
    kinds = {c["name"]: c["kind"] for c in payload["columns"]}
    assert kinds["run_key"] == "string" and kinds["seed"] == "integer" and kinds["directed_pr_auc"] == "number"
    assert len(payload["data"]["run_key"]) == 2
    # Per fit (row index) and bin of width 0.1: the 1->3 edge at 0.83 is true, in bin 8; 2->3 and 1->2 sit in bin 0,
    # one of them true.
    calibration = payload["calibration"]
    assert calibration["bins"] == 10
    by_bin = {}
    for row, b, n, n_true, score_sum in zip(*(calibration[k] for k in ("row", "bin", "n", "n_true", "score_sum"))):
        by_bin.setdefault(row, {})[b] = (n, n_true)
        assert 0 <= score_sum <= n
    assert by_bin[0] == by_bin[1] == {0: (2, 1), 8: (1, 1)}


def test_load_without_an_edges_file_has_no_calibration(rows, tmp_path):
    out = tmp_path / "cohort.parquet"
    cohort.concatenate([rows], out)
    cohort.edges_path(out).unlink()
    assert cohort.load(out)["calibration"] is None
