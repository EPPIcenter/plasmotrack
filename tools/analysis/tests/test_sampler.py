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
def client(golden_analysis, tmp_path):
    from starlette.testclient import TestClient

    from pt_analysis.server import create_app

    return TestClient(create_app(golden_analysis, static_dir=tmp_path / "dist"))


def test_the_endpoint_serves_one_chains_panel(client):
    got = client.get("/api/sampler", params={"chain": "chain2"}).json()
    assert got["chain"] == "chain2"
    assert got["windows"] == [
        {"window": 0, "phase": "burn_in", "first_iteration": 0, "last_iteration": 20},
        {"window": 1, "phase": "sampling", "first_iteration": 21, "last_iteration": 45},
        {"window": 2, "phase": "sampling", "first_iteration": 46, "last_iteration": 70}]
    acceptance = got["acceptance"]
    assert len(acceptance["type"]) == 3 * len(set(acceptance["type"]))  # a row per sampler type per Window
    assert sum(acceptance["targets"]) == 108  # 36 samplers in each of 3 Windows
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


# Per-type aggregation: the panel's default payload is one row per sampler type per Window, columnar, so a Run with
# tens of thousands of samplers stays small. Expected values below are worked by hand from `ACCEPTANCE`.

ACCEPTANCE = [  # (window, sampler, accepted, proposed)
    (0, "Relink/a", 2, 10), (1, "Relink/a", 9, 10),   # rates 0.2, 0.9: mean 0.55
    (0, "Relink/b", 5, 10), (1, "Relink/b", 0, 0),    # rates 0.5, none: mean 0.5
    (0, "Relink/c", 8, 20), (1, "Relink/c", 1, 10),   # rates 0.4, 0.1: mean 0.25
    (0, "Walk/x", 3, 4), (1, "Walk/x", 1, 4),
]


def acceptance_store(path, rows=ACCEPTANCE, variances=None):
    """A store with Windows 0 (burn-in) and 1 (sampling), an `acceptance` table of `rows`, `structure/samplers` for
    them and, with `variances`, a `burn_in_ended` event carrying them."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    windows = [(0, 0, 20, "burn_in"), (1, 21, 45, "sampling")]
    store = write_store(path, windows, [event(20, "burn_in", 0, "burn_in_ended",
                                              {"proposal_variances": variances or {}, "temperatures": [1.0]})])
    phases = {n: phase for n, _, _, phase in windows}
    (path / "acceptance").mkdir()
    for number in phases:
        pq.write_table(pa.Table.from_pylist(
            [{"window": w, "phase": phases[w], "sampler": s, "accepted": a, "proposed": p}
             for w, s, a, p in rows if w == number]), path / "acceptance" / f"window-{number:06d}.parquet")
    (path / "structure").mkdir()
    names = sorted({s for _, s, _, _ in rows})
    pq.write_table(pa.table({"sampler": names, "type": [n.split("/")[0] for n in names],
                             "target": [n.split("/", 1)[1] for n in names]}), path / "structure" / "samplers.parquet")
    return store


def test_acceptance_by_type_sums_each_types_counts_per_window_with_the_spread_of_target_rates(tmp_path):
    got = sampler.acceptance_by_type(acceptance_store(tmp_path / "run"))
    assert got["type"] == ["Relink", "Relink", "Walk", "Walk"]
    assert got["window"] == [0, 1, 0, 1]
    assert got["targets"] == [3, 3, 1, 1]
    assert got["accepted"] == [15, 10, 3, 1]
    assert got["proposed"] == [40, 20, 4, 4]
    assert got["rate"] == pytest.approx([0.375, 0.5, 0.75, 0.25])
    assert got["rate_min"] == pytest.approx([0.2, 0.1, 0.75, 0.25])
    assert got["rate_median"] == pytest.approx([0.4, 0.5, 0.75, 0.25])  # Window 1: b proposed nothing, has no rate
    assert got["rate_max"] == pytest.approx([0.5, 0.9, 0.75, 0.25])


def test_a_types_targets_come_lowest_mean_rate_first_capped_with_the_number_cut(tmp_path):
    got = sampler.acceptance_targets(acceptance_store(tmp_path / "run"), "Relink", limit=2)
    assert (got["type"], got["targets"], got["shown"]) == ("Relink", 3, 2)
    rows = got["acceptance"]
    assert rows["target"] == ["c", "c", "b", "b"]  # mean rates: c 0.25, b 0.5; a (0.55) is cut
    assert rows["sampler"] == ["Relink/c", "Relink/c", "Relink/b", "Relink/b"]
    assert rows["window"] == [0, 1, 0, 1]
    assert rows["accepted"] == [8, 1, 5, 0]
    assert rows["proposed"] == [20, 10, 10, 0]
    assert rows["rate"] == pytest.approx([0.4, 0.1, 0.5, None])


def test_a_types_targets_default_to_a_cap_of_500(tmp_path):
    rows = [(w, f"Relink/{i}", i % 7, 10) for i in range(600) for w in (0, 1)]
    got = sampler.acceptance_targets(acceptance_store(tmp_path / "run", rows), "Relink")
    assert (got["targets"], got["shown"], len(got["acceptance"]["target"])) == (600, 500, 1000)
    assert set(got["acceptance"]["accepted"]) == {0, 1, 2, 3, 4, 5}  # every target accepting 6 of 10 is cut


def test_a_type_with_no_samplers_has_no_targets(tmp_path):
    got = sampler.acceptance_targets(acceptance_store(tmp_path / "run"), "Nope")
    assert (got["targets"], got["shown"], got["acceptance"]["target"]) == (0, 0, [])


VARIANCES = {"Walk/x": [1.0], "Walk/y": [0.5, 2.0], "Relink/b": [3.0], "Relink/c": [4.0, 0.25, 0.5]}


def test_the_panel_summarises_burn_in_end_proposal_variances_by_type(tmp_path):
    end = sampler.panel(acceptance_store(tmp_path / "run", variances=VARIANCES))["burn_in_ended"]
    assert end["iteration"] == 20 and end["temperatures"] == [1.0]
    assert end["proposal_variances"] == {"type": ["Relink", "Walk"], "samplers": [2, 2], "scales": [4, 3],
                                         "min": [0.25, 0.5], "median": [1.75, 1.0], "max": [4.0, 2.0]}


def test_the_panel_aggregates_acceptance_by_type(tmp_path):
    panel = sampler.panel(acceptance_store(tmp_path / "run"))
    assert panel["acceptance"]["type"] == ["Relink", "Relink", "Walk", "Walk"]
    assert panel["burn_in_ended"]["proposal_variances"]["type"] == []


def test_a_types_panel_gives_its_shown_targets_proposal_variances_in_their_order(tmp_path):
    got = sampler.type_panel(acceptance_store(tmp_path / "run", variances=VARIANCES), "Relink", limit=2)
    assert got["acceptance"]["target"] == ["c", "c", "b", "b"]
    assert got["proposal_variances"] == {"target": ["c", "b"], "sampler": ["Relink/c", "Relink/b"],
                                         "variances": [[4.0, 0.25, 0.5], [3.0]]}


def test_the_endpoint_sums_a_types_targets_per_window(client, golden_analysis):
    got = client.get("/api/sampler", params={"chain": "chain2"}).json()["acceptance"]
    rows = sampler.acceptance(golden_analysis.runs["chain2"])
    walk = [i for i, t in enumerate(got["type"]) if t == "ConstrainedContinuousRandomWalk"]
    for i in walk:
        mine = [r for r in rows if r["type"] == got["type"][i] and r["window"] == got["window"][i]]
        assert got["accepted"][i] == sum(r["accepted"] for r in mine)
        assert got["proposed"][i] == sum(r["proposed"] for r in mine)
        assert got["targets"][i] == len(mine)


def test_the_endpoint_serves_one_types_targets_capped(client):
    got = client.get("/api/sampler", params={"chain": "chain2", "type": "ConstrainedContinuousRandomWalk",
                                             "limit": "2"}).json()
    assert got["chain"] == "chain2" and got["type"] == "ConstrainedContinuousRandomWalk"
    assert got["shown"] == 2 and got["targets"] == 12  # 12 random walks in the golden Runs
    assert len(got["acceptance"]["target"]) == 6  # 2 targets × 3 Windows
    assert got["proposal_variances"]["target"] == list(dict.fromkeys(got["acceptance"]["target"]))


def test_the_endpoint_refuses_an_unknown_type_and_a_bad_limit(client):
    response = client.get("/api/sampler", params={"type": "NoSuchSampler"})
    assert response.status_code == 404 and "NoSuchSampler" in response.json()["error"]
    assert client.get("/api/sampler", params={"type": "SALTSampler", "limit": "x"}).status_code == 400
    assert client.get("/api/sampler", params={"type": "SALTSampler", "limit": "0"}).status_code == 400


def test_the_panels_burn_in_ended_event_leaves_out_the_per_sampler_variances(tmp_path):
    events = sampler.panel(acceptance_store(tmp_path / "run", variances=VARIANCES))["events"]
    assert [(e["type"], e["payload"]) for e in events] == [("burn_in_ended", {"temperatures": [1.0]})]
