"""The Store reads a Run's store under the closed-Window rule, against the golden stores."""

from conftest import GOLDEN
from pt_analysis.store import Store, Window


def test_closed_windows_come_from_the_markers(golden_store):
    assert Store(golden_store).closed_windows() == [
        Window(number=0, first_iteration=0, last_iteration=20, phase="burn_in"),
        Window(number=1, first_iteration=21, last_iteration=45, phase="sampling"),
        Window(number=2, first_iteration=46, last_iteration=70, phase="sampling"),
    ]


def test_a_trace_is_read_over_the_sampling_phase_by_default(golden_store):
    rows = Store(golden_store).table("run_scalars").fetchall()
    assert len(rows) == 10


def test_a_trace_can_be_read_over_both_phases(golden_store):
    store = Store(golden_store)
    assert len(store.table("run_scalars", phases=("burn_in", "sampling")).fetchall()) == 15
    assert len(store.table("run_scalars", phases=("burn_in",)).fetchall()) == 5


def test_parent_sets_rows_per_phase():
    single, re = Store(GOLDEN / "single-chain"), Store(GOLDEN / "replica-exchange")
    assert len(single.table("parent_sets").fetchall()) == 90
    assert len(re.table("parent_sets").fetchall()) == 86
    assert len(single.table("parent_sets", phases=("burn_in",)).fetchall()) == 37
    assert len(re.table("parent_sets", phases=("burn_in",)).fetchall()) == 43


def test_parent_set_members_are_read_over_both_phases_by_default(golden_store):
    members = Store(golden_store).table("parent_set_members")
    rows = members.order("infection, set_id, member").fetchall()
    assert rows == [
        ("1", 0, "1_copy"),
        ("2", 0, "2_copy"), ("2", 1, "1"), ("2", 1, "2_copy"),
        ("3", 0, "3_copy"), ("3", 1, "1"), ("3", 1, "3_copy"), ("3", 2, "2"), ("3", 2, "3_copy"),
        ("3", 3, "1"), ("3", 3, "2"), ("3", 3, "3_copy"), ("3", 4, "1"), ("3", 5, "1"), ("3", 5, "2"),
    ]


def test_a_table_with_no_rows_in_the_chosen_windows_is_none(golden_store):
    store = Store(golden_store)
    assert store.table("parent_set_members", phases=("sampling",)) is None
    assert store.table("latent_parent_genotypes") is None


def test_structure_tables_and_manifest(golden_store):
    store = Store(golden_store)
    assert store.structure("infections").order("infection").fetchall() == [
        ("1", "1_copy", True), ("2", "2_copy", True), ("3", "3_copy", False)]
    assert len(store.structure("loci").fetchall()) == 20
    assert len(store.structure("samplers").fetchall()) == 36
    assert store.manifest["seed"] == 42
    assert store.manifest["window"] == 25


def test_rows_are_dicts_in_iteration_order(golden_store):
    rows = Store(golden_store).rows("run_scalars")
    assert [r["iteration"] for r in rows] == [25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
    assert all(r["phase"] == "sampling" for r in rows)
    assert rows[0]["mean_coi"] > 0
    assert Store(golden_store).rows("latent_parent_genotypes") == []


def test_draws_count_the_written_sampling_iterations(golden_store):
    assert Store(golden_store).draws() == 10


def test_sampling_time_sums_the_window_closed_events():
    assert Store(GOLDEN / "single-chain").sampling_time() == {"sampling_seconds": 0.058744, "samples": 50}
    assert Store(GOLDEN / "replica-exchange").sampling_time() == {"sampling_seconds": 0.260588, "samples": 50}
