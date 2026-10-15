"""A Window exists only if its marker exists: files of an unmarked Window are ignored by every reader entry point."""

import shutil

import duckdb
import pytest
from conftest import GOLDEN

from pt_analysis.store import Store


@pytest.fixture
def store_with_unmarked_window(tmp_path):
    """The single-chain golden store plus `window-000003.parquet` files without a marker, as a Run killed
    mid-commit leaves them, and an in-flight file."""
    path = tmp_path / "store"
    shutil.copytree(GOLDEN / "single-chain", path)
    for table in ["run_scalars", "parent_sets", "events", "infections"]:
        src, dst = path / table / "window-000002.parquet", path / table / "window-000003.parquet"
        duckdb.sql(f"""COPY (SELECT * REPLACE (3 AS "window", iteration + 25 AS iteration) FROM '{src}') TO '{dst}'""")
    duckdb.sql(f"COPY (SELECT '3' AS infection, 6 AS set_id, '9' AS member) "
               f"TO '{path / 'parent_set_members' / 'window-000003.parquet'}'")
    shutil.copy(path / "run_scalars" / "window-000003.parquet", path / "run_scalars" / ".in-flight.parquet")
    # The fixture has teeth: a reader that globbed the directory would see Window 3.
    seen = duckdb.sql(f"""SELECT DISTINCT "window" FROM '{path / 'run_scalars' / 'window-*.parquet'}'""").fetchall()
    assert (3,) in seen
    return Store(path)


def test_closed_windows_ignore_the_unmarked_window(store_with_unmarked_window):
    assert [w.number for w in store_with_unmarked_window.closed_windows()] == [0, 1, 2]


def test_tables_ignore_the_unmarked_window(store_with_unmarked_window):
    store = store_with_unmarked_window
    for table in ["run_scalars", "parent_sets", "events", "infections"]:
        windows = store.table(table, phases=("burn_in", "sampling")).aggregate('list(DISTINCT "window")').fetchone()[0]
        assert sorted(windows) == [0, 1, 2], table
    assert len(store.table("run_scalars").fetchall()) == 10
    assert len(store.table("parent_set_members").fetchall()) == 15


def test_rows_draws_and_sampling_time_ignore_the_unmarked_window(store_with_unmarked_window):
    store = store_with_unmarked_window
    assert [r["iteration"] for r in store.rows("run_scalars")][-1] == 70
    assert store.draws() == 10
    assert store.sampling_time() == Store(GOLDEN / "single-chain").sampling_time()
