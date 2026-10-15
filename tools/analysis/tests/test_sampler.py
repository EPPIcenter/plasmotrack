"""The sampler panel's data (U4): acceptance per sampler and swap acceptance per temperature pair per Window, the
temperature ladder over adaptation, the burn-in-end proposal variances and ladder, and the Run events marked on the
panel's time axes. Golden-store values are in the exploration notes and `doc/store-layout.md`."""

import pytest
from conftest import GOLDEN

from pt_analysis import sampler
from pt_analysis.store import Store


def test_acceptance_has_a_row_per_sampler_per_window_over_both_phases(store):
    rows = sampler.acceptance(store)
    assert len(rows) == 108
    assert {r["window"] for r in rows} == {0, 1, 2}
    assert len({r["sampler"] for r in rows}) == 36
    assert {(r["window"], r["phase"]) for r in rows} == {(0, "burn_in"), (1, "sampling"), (2, "sampling")}


def test_acceptance_rows_carry_type_target_window_range_and_rate(golden_store):
    rows = sampler.acceptance(Store(golden_store))
    row = next(r for r in rows if r["window"] == 0 and r["sampler"] == "ConstrainedContinuousRandomWalk/mean_coi")
    assert row["type"] == "ConstrainedContinuousRandomWalk"
    assert row["target"] == "mean_coi"
    assert (row["first_iteration"], row["last_iteration"]) == (0, 20)
    assert row["proposed"] == 20
    assert row["rate"] == pytest.approx(row["accepted"] / 20)
    if golden_store.name == "single-chain":
        assert row["accepted"] == 14  # exploration notes: `0 burn_in …/mean_coi 14 20`


def test_a_sampler_that_proposed_nothing_has_no_rate(store):
    rows = sampler.acceptance(store)
    row = next(r for r in rows if r["window"] == 0 and r["sampler"] == "CloneOrderSwapSampler/ordering/1")
    assert (row["accepted"], row["proposed"], row["rate"]) == (0, 0, None)


def test_acceptance_rows_are_ordered_by_type_target_and_window(store):
    rows = sampler.acceptance(store)
    keys = [(r["type"], r["target"], r["window"]) for r in rows]
    assert keys == sorted(keys)


def replica_exchange():
    return Store(GOLDEN / "replica-exchange")


def single_chain():
    return Store(GOLDEN / "single-chain")


def test_swaps_have_a_row_per_temperature_pair_per_window():
    rows = sampler.swaps(replica_exchange())
    assert len(rows) == 9  # 3 pairs of 4 temperature slots × 3 Windows
    assert rows[3] == {"window": 1, "phase": "sampling", "first_iteration": 21, "last_iteration": 45,
                       "temperature_pair": 0, "accepted": 12, "attempted": 13, "rate": pytest.approx(12 / 13)}
    assert [(r["window"], r["temperature_pair"]) for r in rows] == [(w, k) for w in range(3) for k in range(3)]


def test_a_run_with_one_chain_has_no_swaps():
    assert sampler.swaps(single_chain()) == []


def test_the_ladder_is_read_from_temperatures_adapted_events():
    ladder = sampler.ladder(replica_exchange())
    assert len(ladder) == 1
    step = ladder[0]
    assert (step["iteration"], step["phase"], step["applied"]) == (16, "burn_in", True)
    assert step["temperatures"] == pytest.approx([1, 0.8486, 0.7431, 0.625], abs=1e-4)
    assert step["barriers"] == pytest.approx([0.11527, 0.17173, 0.15777], abs=1e-5)


def test_a_run_with_one_chain_has_no_ladder():
    assert sampler.ladder(single_chain()) == []


@pytest.mark.parametrize("name, temperatures", [("single-chain", [1.0]),
                                                ("replica-exchange", [1, 0.8486, 0.7431, 0.625])])
def test_burn_in_end_gives_the_final_proposal_variances_and_ladder(name, temperatures):
    end = sampler.burn_in_end(Store(GOLDEN / name))
    assert end["iteration"] == 20
    assert end["temperatures"] == pytest.approx(temperatures, abs=1e-4)
    variances = end["proposal_variances"]
    assert len(variances) == 14  # samplers with a proposal scale: 12 random walks and 2 SALT samplers
    assert variances["ConstrainedContinuousRandomWalk/mean_coi"] == [1.0]
    assert variances["SALTSampler/allele_frequencies/AS1"] == [1.0] * 8  # one per Allele
    assert len(variances["SALTSampler/allele_frequencies/AS11"]) == 12


def test_the_golden_runs_mark_burn_in_end_and_run_end(golden_store):
    events = sampler.marked_events(Store(golden_store))
    assert [(e["iteration"], e["type"]) for e in events] == [(20, "burn_in_ended"), (70, "run_ended")]
    assert events[1]["payload"]["reason"] == "completed"
    assert events[1]["window"] == 2 and events[1]["phase"] == "sampling"


def write_store(path, windows, events):
    """A store with closed Windows `windows` ([(number, first, last, phase)]) and an `events` table, one file per
    Window; event rows name their Window."""
    import json

    import pyarrow as pa
    import pyarrow.parquet as pq

    (path / "_closed").mkdir(parents=True)
    (path / "events").mkdir()
    for number, first, last, phase in windows:
        (path / "_closed" / f"window-{number:06d}").write_text(json.dumps(
            {"first_iteration": first, "last_iteration": last, "phase": phase, "window": number}))
    for number in {e["window"] for e in events}:
        rows = [{**e, "payload": json.dumps(e["payload"])} for e in events if e["window"] == number]
        pq.write_table(pa.Table.from_pylist(rows), path / "events" / f"window-{number:06d}.parquet")
    return Store(path)


def event(iteration, phase, window, type, payload):
    return {"iteration": iteration, "phase": phase, "window": window, "wall_time": "2026-10-14T12:00:00.000Z",
            "type": type, "payload": payload}


def test_non_finite_values_and_interrupts_are_marked_with_their_payloads(tmp_path):
    store = write_store(tmp_path / "run", [(0, 0, 20, "burn_in"), (1, 21, 33, "sampling")], [
        event(0, "burn_in", 0, "run_started", {"burn_in": 20, "samples": 50}),
        event(7, "burn_in", 0, "non_finite_value", {"temperature_slot": 2, "term": "likelihood", "value": "-inf"}),
        event(20, "burn_in", 0, "burn_in_ended", {"proposal_variances": {}, "temperatures": [1.0]}),
        event(20, "burn_in", 0, "window_closed", {"sampling_iterations": 20, "seconds": 1.0}),
        event(33, "sampling", 1, "interrupted", {"signal": "SIGINT", "second": False}),
        event(33, "sampling", 1, "run_ended", {"reason": "interrupted", "seconds": 2.0}),
        # A Window with no marker: a Run killed mid-commit. Its events are not read.
        event(40, "sampling", 2, "non_finite_value", {"temperature_slot": 0, "term": "prior", "value": "nan"}),
    ])
    events = sampler.marked_events(store)
    assert [(e["iteration"], e["type"]) for e in events] == [
        (7, "non_finite_value"), (20, "burn_in_ended"), (33, "interrupted"), (33, "run_ended")]
    assert events[0] == {"iteration": 7, "phase": "burn_in", "window": 0, "wall_time": "2026-10-14T12:00:00.000Z",
                         "type": "non_finite_value",
                         "payload": {"temperature_slot": 2, "term": "likelihood", "value": "-inf"}}
    assert events[2]["payload"] == {"signal": "SIGINT", "second": False}


def test_a_resumed_run_marks_where_it_resumed(tmp_path):
    store = write_store(tmp_path / "run", [(0, 0, 20, "burn_in"), (1, 21, 30, "sampling"), (2, 31, 45, "sampling")], [
        event(30, "sampling", 1, "run_ended", {"reason": "interrupted", "seconds": 2.0}),
        event(30, "sampling", 2, "run_resumed", {"iteration": 30, "window": 1, "checkpoint": "c", "samples": 25}),
    ])
    assert [e["type"] for e in sampler.marked_events(store)] == ["run_ended", "run_resumed"]
    assert sampler.ladder(store) == [] and sampler.burn_in_end(store) is None


@pytest.fixture
def client(tmp_path):
    import shutil

    from starlette.testclient import TestClient
    from test_analysis import write_manifest

    from pt_analysis.analysis import Analysis
    from pt_analysis.server import create_app

    for name in ("single-chain", "replica-exchange"):
        shutil.copytree(GOLDEN / name, tmp_path / name)
    analysis = Analysis.open(write_manifest(tmp_path / "analysis.json",
                                            runs={"chain1": "single-chain", "chain2": "replica-exchange"}))
    return TestClient(create_app(analysis, static_dir=tmp_path / "dist"))


def test_the_endpoint_serves_one_chains_panel(client):
    got = client.get("/api/sampler", params={"chain": "chain2"}).json()
    assert got["chain"] == "chain2"
    assert got["windows"] == [
        {"window": 0, "phase": "burn_in", "first_iteration": 0, "last_iteration": 20},
        {"window": 1, "phase": "sampling", "first_iteration": 21, "last_iteration": 45},
        {"window": 2, "phase": "sampling", "first_iteration": 46, "last_iteration": 70}]
    assert len(got["acceptance"]) == 108
    assert got["swaps"][3]["accepted"] == 12 and got["swaps"][3]["attempted"] == 13
    assert got["ladder"][0]["iteration"] == 16
    assert got["burn_in_ended"]["temperatures"] == pytest.approx([1, 0.8486, 0.7431, 0.625], abs=1e-4)
    assert [e["type"] for e in got["events"]] == ["burn_in_ended", "run_ended"]


def test_the_endpoint_defaults_to_the_first_chain_and_refuses_an_unknown_one(client):
    got = client.get("/api/sampler").json()
    assert got["chain"] == "chain1" and got["swaps"] == [] and got["ladder"] == []
    response = client.get("/api/sampler", params={"chain": "chain9"})
    assert response.status_code == 404
    assert "chain9" in response.json()["error"]
